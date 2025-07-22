#ifndef GEL_HMesh_RsR2_hpp
#define GEL_HMesh_RsR2_hpp
#pragma once

// TODO: unnecessary imports will negatively affect compile times
#include <random>
#include <vector>
#include <ranges>

#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/etf.h>
#include <GEL/Geometry/KDTree.h>

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

using Tree = Geometry::KDTree<Point, NodeID>;
using Record = Geometry::KDTreeRecord<Point, NodeID>;

struct NeighborInfo {
    NodeID id;
    double distance; // the added precision doesn't justify the usage of doubles

    static constexpr NodeID invalid_id = -1;
    NeighborInfo() = delete;

    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
};

using NeighborArray = std::vector<NeighborInfo>;
using NeighborMap = std::vector<NeighborArray>;

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

/// A trivial wrapper over bool to avoid the std::vector<bool> specialization
struct Boolean {
    bool inner;
};

/*Graph definition. The RsR graph here is integrated with the rotation system based on AMGraph*/
struct Vertex {
    NodeID id = 0;
    int normal_rep = -1;
    Vec3 coords = Vec3(0., 0., 0.);
    Vec3 normal = Vec3(0., 0., 0.);

    bool operator==(const Vertex& rhs) const { return id == rhs.id; }

    struct Neighbor {
        double angle;
        uint v;
        mutable uint tree_id = 0;

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
    double total_edge_length = 0.;
    int exp_genus = -1;
    // These should be function parameters, not fields
    bool is_euclidean = false;
    bool is_final = false;
    size_t current_no_edges = 0;

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

        current_no_edges++;
        m_edges[id].weight = weight;
        m_edges[id].source = source;
        m_edges[id].target = target;
        this->total_edge_length += weight;
        insert_neighbor(source, target);
        insert_neighbor(target, source);

        return id;
    }

    NodeID add_node(const Vec3& p, const Vec3& in_normal = Vec3(0., 0., 0.))
    {
        const NodeID n = AMGraph::add_node();
        m_vertices[n] = Vertex{.id = n, .normal_rep = -1, .coords = p, .normal = in_normal};
        return n;
    }

    void init(const std::vector<Point>& vertices, const std::vector<Vec3>& normals)
    {
        for (size_t i = 0; i < vertices.size(); i++) {
            const NodeID id = this->add_node(vertices[i]);
            m_vertices[id].normal = normals[i];
        }
    }

    void init(const std::vector<Point>& vertices)
    {
        for (const auto& vertex : vertices) {
            this->add_node(vertex);
        }
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

auto collapse_points(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    size_t max_iterations
) -> Collapse;
} // namespace HMesh::RSR


#endif // GEL_HMesh_RsR2_hpp
