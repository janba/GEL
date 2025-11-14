#ifndef GEL_RSR_H
#define GEL_RSR_H
#pragma once

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/obj_load.h>
#include <GEL/Geometry/etf.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/normal.h>

namespace HMesh
{
namespace detail
{
    using namespace CGLA;
    using namespace Geometry;
    using NodeID = AMGraph::NodeID;

    typedef CGLA::Vec3d Vector;
    using Point = Vector;
    typedef std::pair<NodeID, NodeID> m_Edge;

    double cal_radians_3d(const Vector& branch, const Vector& normal);

    double cal_radians_3d(const Vector& branch_vec, const Vector& normal,
                          const Vector& ref_vec);

    /*Graph definition. The RsR graph here is integrated with rotation system based on AMGraph*/

    struct Vertex {
        int id = 0;
        int normal_rep = -1;
        bool colored = false;
        Vector coords = Vector(0., 0., 0.);
        Vector normal = Vector(0., 0., 0.);
        std::vector<bool> faceExist;
        float distance;

        struct Neighbor {
            double angle;
            uint v;
            mutable uint tree_id = 0;
            mutable bool faceExist = false;

            Neighbor(const Vertex& u, const Vertex& v, uint id)
            {
                this->v = id;
                //std::cout << v.coords << std::endl;
                //std::cout << u.coords << std::endl;
                //std::cout << cal_radians_3d(v.coords - u.coords, u.normal) << std::endl;
                this->angle = cal_radians_3d(v.coords - u.coords, u.normal);
            }

            bool operator <(const Neighbor& rhs) const
            {
                return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
            }

            bool operator<(Neighbor& rhs) const
            {
                return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
            }
        };

        std::set<Neighbor> ordered_neighbors;
    };

    struct Edge {
        NodeID source = -1;
        NodeID target = -1;
        double weight = 0.;
        int ref_time = 0;
        int count_weight = 1;
    };

    typedef Vertex::Neighbor Neighbor;

    class SimpGraph : public AMGraph {
    public:
        Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

        EdgeID connect_nodes(NodeID source, NodeID target, float weight = 0.)
        {
            // std::cout << __FILE__ << " " << __LINE__ << " " << source << " " << target << std::endl;
            EdgeID id = AMGraph::connect_nodes(source, target);
            // std::cout << __FILE__ << " " << __LINE__ << " " << id << std::endl;
            m_edges[id].weight = weight;
            return id;
        }

        double get_weight(NodeID n1, NodeID n2) const
        {
            double output = m_edges[find_edge(n1, n2)].weight;
            return output;
        }

        /** Disconnect nodes. This operation removes the edge from the edge maps of the two formerly connected
             vertices, but the number of edges reported by the super class AMGraph is not decremented, so the edge is only
             invalidated. Call cleanup to finalize removal. */
        void disconnect_nodes(NodeID n0, NodeID n1)
        {
            if (valid_node_id(n0) && valid_node_id(n1)) {
                edge_map[n0].erase(n1);
                edge_map[n1].erase(n0);
            }
        }
    };

    class RSGraph : public AMGraph {
    public:
        double total_edge_length = 0.;
        int face_loop_id = 0;
        bool isEuclidean = false;
        bool isFinal = false;
        int exp_genus = -1;
        ETF etf;
        int current_no_edges = 0;

        Util::AttribVec<NodeID, Vertex> m_vertices;
        Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

        /// Compute sqr distance between two nodes - not necessarily connected.
        double dist(NodeID n0, NodeID n1) const
        {
            if (valid_node_id(n0) && valid_node_id(n1))
                return (m_vertices[n0].coords - m_vertices[n1].coords).length();
            else
                return CGLA::CGLA_NAN;
        }

        /// Compute the average edge length
        double cal_average_edge_length()
        {
            double sum_len = this->total_edge_length;
            return sum_len / current_no_edges;
        }

        void remove_edge(NodeID source, NodeID target)
        {
            this->total_edge_length -= this->m_edges[edge_map[source][target]].weight;
            if (valid_node_id(source) && valid_node_id(target)) {
                edge_map[source].erase(target);
                edge_map[target].erase(source);
            }
            current_no_edges--;
            return;
        }

        void remove_neighbor(NodeID root, NodeID neighbor)
        {
            auto& u = m_vertices[root];
            auto& v = m_vertices[neighbor];
            u.ordered_neighbors.erase(Neighbor(u, v, neighbor));
            return;
        }

        void insert_neighbor(NodeID root, NodeID neighbor)
        {
            const auto& u = m_vertices[root];
            const auto& v = m_vertices[neighbor];
            m_vertices[root].ordered_neighbors.insert(Neighbor(u, v, neighbor));
            //std::cout << Neighbor(u, v, neighbor).angle << std::endl;
            return;
        }

        EdgeID add_edge(NodeID source, NodeID target, float weight = 0.)
        {
            EdgeID id = this->connect_nodes(source, target);
            if (id != InvalidEdgeID) {
                current_no_edges++;
                m_edges[id].weight = weight;
                m_edges[id].source = source;
                m_edges[id].target = target;
                this->total_edge_length += weight;
                insert_neighbor(source, target);
                insert_neighbor(target, source);
            } else {
                /*std::cout << "weird" << std::endl;
                std::cout << valid_node_id(source) << std::endl;
                std::cout << valid_node_id(target) << std::endl;
                std::cout << source << std::endl;
                std::cout << target << std::endl;
                std::cout << (find_edge(target, source) == InvalidEdgeID)
                    << std::endl;*/
            }

            return id;
        }

        NodeID add_node(const Vector& p)
        {
            NodeID n = AMGraph::add_node();
            Vertex v;
            v.id = n;
            v.coords = p;
            m_vertices[n] = v;
            return n;
        }

        NodeID add_node(const Vector& p, const Vector& in_normal)
        {
            NodeID n = AMGraph::add_node();
            Vertex v;
            v.id = n;
            v.coords = p;
            v.normal = in_normal;
            m_vertices[n] = v;
            return n;
        }

        void init(const std::vector<Point>& vertices, const std::vector<Vector> normals)
        {
            for (int i = 0; i < vertices.size(); i++) {
                NodeID id = this->add_node(vertices[i]);
                m_vertices[id].normal = normals[i];
            }
        }

        void init(const std::vector<Point>& vertices)
        {
            for (int i = 0; i < vertices.size(); i++) {
                this->add_node(vertices[i]);
            }
        }

        void init(int no_vertex)
        {
            for (int i = 0; i < no_vertex; i++) {
                AMGraph::add_node();
            }
        }

        void get_node_set(NodeSet& sets)
        {
            for (NodeID i = 0; i < edge_map.size(); i++) {
                sets.insert(i);
            }
            return;
        }
    };

    typedef Geometry::KDTree<Point, NodeID> Tree;
    typedef Geometry::KDTreeRecord<Point, NodeID> Record;

    void kNN_search(const Point&, const Tree&, int,
                    std::vector<NodeID>&, std::vector<double>&,
                    double last_dist = INFINITY, bool isContain = true);

    void NN_search(const Point&, const Tree&, double,
                   std::vector<NodeID>&, std::vector<double>&, bool isContain = true);

    float find_components(std::vector<Point>&,
                          std::vector<std::vector<Point>>&, std::vector<Point>&,
                          std::vector<std::vector<Point>>&, std::vector<Vector>&,
                          std::vector<std::vector<Vector>>&, const Tree&, float, float);

    void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
                    const std::vector<Vector>& normals, const Tree& kdTree, SimpGraph& dist_graph,
                    std::vector<float>& max_length, std::vector<float>& pre_max_length, float cross_conn_thresh);

    int find_shortest_path(const RSGraph& mst, NodeID start, NodeID target,
                           int threshold, std::vector<NodeID>& path);

    void weighted_smooth(const std::vector<Point>& vertices,
                         std::vector<Point>& smoothed_v, const std::vector<Vector>& normals,
                         const Tree& kdTree, float diagonal_length);

    void estimate_normal(const std::vector<Point>& vertices,
                         const Tree& kdTree, std::vector<Vector>& normals,
                         std::vector<NodeID>& zero_normal_id, float& diagonal_length);

    void minimum_spanning_tree(const SimpGraph& g, NodeID root,
                               RSGraph& gn, std::vector<Vector>& normals, std::vector<Point>& vertices);

    void minimum_spanning_tree(const SimpGraph& g, NodeID root, SimpGraph& gn);

    void correct_normal_orientation(std::vector<Point>& in_smoothed_v,
                                    Tree& kdTree, std::vector<Vector>& normals);

    bool register_face(RSGraph& mst, NodeID v1, NodeID v2, std::vector<std::vector<int>>& faces,
                       Tree& KDTree, float edge_length);

    void add_face(RSGraph& G, std::vector<NodeID>& item,
                  std::vector<std::vector<NodeID>>& faces);

    void connect_handle(const std::vector<Point>& smoothed_v, Tree& KDTree,
                        RSGraph& mst, std::vector<NodeID>& connected_handle_root,
                        std::vector<int>& betti);

    // Timer

    class RsR_Timer {
    public:
        RsR_Timer() {}
        std::vector<long> times;
        std::vector<std::string> descriptions;
        std::vector<std::chrono::high_resolution_clock::time_point> starts;
        std::vector<std::chrono::high_resolution_clock::time_point> ends;

        void create(std::string name)
        {
            times.push_back(0);
            descriptions.push_back(name);
            //idx_map.insert(std::pair<std::string, int>(name, times.size() - 1));
            starts.push_back(std::chrono::high_resolution_clock::now());
            ends.push_back(std::chrono::high_resolution_clock::now());
        }

        void start(std::string name)
        {
            //int idx = idx_map[name];
            int idx = -1;
            for (int i = 0; i < descriptions.size(); i++) {
                if (descriptions[i] == name) {
                    idx = i;
                    break;
                }
            }
            starts[idx] = std::chrono::high_resolution_clock::now();
        }

        void end(std::string name)
        {
            //int idx = idx_map[name];
            int idx = -1;
            for (int i = 0; i < descriptions.size(); i++) {
                if (descriptions[i] == name) {
                    idx = i;
                    break;
                }
            }
            ends[idx] = std::chrono::high_resolution_clock::now();
            times[idx] +=
                std::chrono::duration_cast<std::chrono::seconds>(ends[idx] - starts[idx]).count();
        }

        void show()
        {
            std::cout << "Time Statistics" << std::endl;
            std::cout << std::string(20, '=') << std::endl;
            for (int i = 0; i < times.size(); i++) {
                std::cout << "Spent " << double(times[i])
                    << " seconds on " << descriptions[i] << std::endl;
            }
        }

        long log(std::string name)
        {
            int idx = -1;
            for (int i = 0; i < descriptions.size(); i++) {
                if (descriptions[i] == name) {
                    idx = i;
                    break;
                }
            }
            return times[idx];
        }
    };

    // Face Loop

    void init_face_loop_label(RSGraph& g);

    const Neighbor& successor(const RSGraph& g,
                              const NodeID& root,
                              const NodeID& branch);


    const Neighbor& predecessor(const RSGraph& g,
                                const NodeID& root,
                                const NodeID& branch);

    void maintain_face_loop(RSGraph& g,
                            const NodeID source, const NodeID target);

    const Neighbor& get_neighbor_info(const RSGraph& g,
                                      const NodeID& root, const NodeID& branch);

    // Utils
    void showProgressBar(float progress);

    Vector projected_vector(Vector& input, Vector& normal);

    void find_common_neighbor(NodeID neighbor, NodeID root,
                              std::vector<NodeID>& share_neighbor, RSGraph& g);

    // Algorithm

    bool geometry_check(RSGraph& mst, m_Edge& candidate,
                        Tree& kdTree);

    bool Vanilla_check(RSGraph& mst, m_Edge& candidate,
                       Tree& kdTree);

    bool isIntersecting(RSGraph& mst, NodeID v1,
                        NodeID v2, NodeID v3, NodeID v4);

    bool routine_check(RSGraph& mst, std::vector<NodeID>& triangle);

    void reset_static();
}


/**
 * @brief Reconstructs a single mesh manifold from input points and normals.
 *
 * This function performs surface reconstruction on a set of input points and normals,
 * producing an output HMesh::Manifold. The reconstruction allows control over several algorithmic
 * parameters such as neighborhood size, radius, angle threshold, and sample count. The difference
 * between Euclidean and projected ditance is that the latter is more resilient to noise. For noise-free data,
 * Euclidean distance is preferred.
 *
 * @param[out] output         The resulting reconstructed manifold.
 * @param[in]  org_vertices   The original input vertices (points) to reconstruct from.
 * @param[in]  org_normals    The corresponding normals for each input vertex.
 * @param[in]  in_isEuclidean True means use Euclidean distance rather than projected (default: false).
 * @param[in]  in_genus       Expected genus of the output surface (default: -1 means auto-detect).
 * @param[in]  in_k           Neighborhood size parameter (default: 70).
 * @param[in]  in_r           Radius parameter for local operations (default: 20).
 * @param[in]  in_theta       Angle threshold parameter in degrees (default: 60).
 * @param[in]  in_n           Number of samples or iterations (default: 50).
 *
 * @note
 *   - The function modifies the output manifold in place.
 *   - Input vectors org_vertices and org_normals must be of the same length.
 *   - Algorithmic parameters may need tuning for different datasets.
 */
void reconstruct_single(HMesh::Manifold& output, std::vector<Vec3d>& org_vertices,
                        std::vector<Vec3d>& org_normals, bool in_isEuclidean = false, int in_genus = -1,
                        int in_k = 70, int in_r = 20, int in_theta = 60, int in_n = 50);
}

#endif
