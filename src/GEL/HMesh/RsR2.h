#ifndef GEL_HMesh_RsR2_hpp
#define GEL_HMesh_RsR2_hpp
#pragma once

#include <vector>

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/etf.h>

#include <GEL/HMesh/Collapse.h>

#include <GEL/Util/ParallelAdapters.h>

/// @namespace HMesh::RSR
namespace HMesh::RSR
{
using namespace CGLA;
using namespace Geometry;

inline const auto InvalidNodeID = AMGraph::InvalidNodeID;
inline const auto InvalidEdgeID = AMGraph::InvalidEdgeID;

using NodeID = AMGraph::NodeID;
using FaceType = std::array<NodeID, 3>;

template <typename T>
using OrderedSet = std::set<T>;

template <typename T>
using UnorderedSet = std::unordered_set<T>;

template <typename K, typename V, typename Hash>
using Map = std::unordered_map<K, V, Hash>;

using Vec3 = Vec3d;
using Point = Vec3;
using TEdge = std::pair<NodeID, NodeID>;

double cal_radians_3d(const Vec3& branch, const Vec3& normal);

double cal_radians_3d(const Vec3& branch_vec, const Vec3& normal,
                      const Vec3& ref_vec);

enum struct Distance {
    EUCLIDEAN,
    NEIGHBORS
};

///
/// TODO: documentation
struct RsROpts {
    int32_t genus = -1;
    int32_t k = 70;
    double r = 20;
    double theta = 60;
    int32_t n = 50;
    Distance dist = Distance::NEIGHBORS;

    bool is_face_normal = true;
    bool is_face_loop = true;
};

/*Graph definition. The RsR graph here is integrated with the rotation system based on AMGraph*/
struct Vertex {
    NodeID id = 0;
    using NormalRep = std::int64_t;
    static constexpr NormalRep InvalidNormalRep = -1;
    static constexpr NormalRep CollisionRep = -2;
    NormalRep normal_rep = InvalidNormalRep;
    Vec3 coords = Vec3(0., 0., 0.);
    Vec3 normal = Vec3(0., 0., 0.);

    bool operator==(const Vertex& rhs) const { return id == rhs.id; }

    struct Neighbor {
        double angle;
        uint v;
        uint tree_id = 0;

        Neighbor(const Vertex& u, const Vertex& v, const uint id)
        {
            this->v = id;
            this->angle = cal_radians_3d(v.coords - u.coords, u.normal);
        }

        friend size_t hash_value(const Neighbor& p)
        {
            return std::hash<uint>()(p.v);
        }

        std::weak_ordering operator<=>(const Neighbor& rhs) const
        {
            if (this->v == rhs.v) return std::weak_ordering::equivalent;
            if (this->angle < rhs.angle) return std::weak_ordering::less;
            else return std::weak_ordering::greater;
        }
    };

    // TODO: mutable elements inside std::set is potentially unsound
    OrderedSet<Neighbor> ordered_neighbors;
};


struct Edge {
    NodeID source = InvalidNodeID;
    NodeID target = InvalidNodeID;
    double weight = 0.;
    // ?
    int ref_time = 0;
};

using Neighbor = Vertex::Neighbor;

class SimpGraph : public AMGraph {
public:
    ::Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

    EdgeID connect_nodes(const NodeID source, const NodeID target, const double weight = 0.)
    {
        const EdgeID id = AMGraph::connect_nodes(source, target);
        m_edges[id].weight = weight;
        return id;
    }

    [[nodiscard]] double get_weight(const NodeID n1, const NodeID n2) const
    {
        return m_edges[find_edge(n1, n2)].weight;
    }

    /** Disconnect nodes. This operation removes the edge from the edge maps of the two formerly connected
         vertices, but the number of edges reported by the super class AMGraph is not decremented, so the edge is only
         invalidated. Call cleanup to finalize removal. */
    void disconnect_nodes(const NodeID n0, const NodeID n1)
    {
        if (valid_node_id(n0) && valid_node_id(n1)) {
            edge_map[n0].erase(n1);
            edge_map[n1].erase(n0);
        }
    }
};

class RSGraph : public AMGraph {
public:
    ETF etf;
    int exp_genus = -1;

    ::Util::AttribVec<NodeID, Vertex> m_vertices;
    ::Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

    void insert_neighbor(const NodeID root, const NodeID neighbor)
    {
        const auto& u = m_vertices[root];
        const auto& v = m_vertices[neighbor];
        m_vertices[root].ordered_neighbors.insert(Neighbor(u, v, neighbor));
    }

    EdgeID add_edge(const NodeID source, const NodeID target, const double weight = 0.)
    {
        const EdgeID id = this->connect_nodes(source, target);
        GEL_ASSERT_NEQ(id, InvalidEdgeID);

        m_edges[id].weight = weight;
        m_edges[id].source = source;
        m_edges[id].target = target;
        insert_neighbor(source, target);
        insert_neighbor(target, source);

        return id;
    }

    NodeID add_node(const Vec3& p, const Vec3& in_normal = Vec3(0., 0., 0.))
    {
        const NodeID n = AMGraph::add_node();
        m_vertices[n] = Vertex{.id = n, .normal_rep = Vertex::InvalidNormalRep, .coords = p, .normal = in_normal};
        return n;
    }

    /// @brief Get last neighbor information
    /// @param root: root vertex index
    /// @param branch: current outgoing branch
    ///
    /// @return reference to last neighbor struct
    [[nodiscard]]
    const Neighbor& predecessor(const NodeID& root, const NodeID& branch) const
    {
        const auto& u = m_vertices.at(root);
        const auto& v = m_vertices.at(branch);
        auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
        // TODO: Issue #105
        // GEL_ASSERT_NEQ(iter, u.ordered_neighbors.end()); // no predecessor exists
        if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
        return (*(std::prev(iter)));
    }

    Neighbor& predecessor(const NodeID& root, const NodeID& branch)
    {
        const auto& u = m_vertices.at(root);
        const auto& v = m_vertices.at(branch);
        auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
        // TODO: Issue #105
        //GEL_ASSERT_NEQ(iter, u.ordered_neighbors.end()); // no predecessor exists
        if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
        return const_cast<Neighbor&>(*(std::prev(iter)));
    }

    /// @brief Get the next neighbor information
    ///
    /// @param root: root vertex index
    /// @param branch: current outgoing branch
    ///
    /// @return reference to the next neighbor struct
    [[nodiscard]]
    const Neighbor& successor(const NodeID& root, const NodeID& branch) const
    {
        const auto& u = m_vertices.at(root);
        const auto& v = m_vertices.at(branch);
        auto iter = u.ordered_neighbors.upper_bound(Neighbor(u, v, branch));
        // TODO: Issue #105
        //GEL_ASSERT_NEQ(iter, u.ordered_neighbors.end()); // no successor exists
        if (iter == u.ordered_neighbors.end()) iter = u.ordered_neighbors.begin(); // Wrap around
        return (*iter); // This is honestly not good practice - ONLY modification of tree_id
    }
};

// Algorithm

/// Convert point cloud to a Manifold
/// @param vertices vertices of the point cloud
/// @param normals normals of the point cloud or empty vector
/// @param opts collapse options
/// @return reconstructed manifold mesh
auto point_cloud_to_mesh(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
                         const RsROpts& opts) -> ::HMesh::Manifold;

struct NormalEstimationResult {
    std::vector<Point> vertices;
    std::vector<Vec3> normals;
    std::vector<Vec3> smoothed_v;
};

auto point_cloud_normal_estimate(const std::vector<Point>& vertices,
                                 const std::vector<Vec3>& normals, bool isEuclidean) -> NormalEstimationResult;


auto point_cloud_collapse_reexpand(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const RsROpts& opts,
    int max_iterations,
    bool reexpand) -> ::HMesh::Manifold;
} // namespace HMesh::RSR


#endif // GEL_HMesh_RsR2_hpp
