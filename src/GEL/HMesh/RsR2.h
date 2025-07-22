#ifndef GEL_HMesh_RsR2_hpp
#define GEL_HMesh_RsR2_hpp
#pragma once

#include <vector>

#include <GEL/Util/AssociativeContainers.h>

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/etf.h>

#include <GEL/HMesh/Collapse.h>

/// @brief Rotation System Reconstruction
namespace HMesh::RSR
{
using Vec3 = CGLA::Vec3d;
using Point = Vec3;
using NodeID = AMGraph::NodeID;

using uint = Geometry::uint;

inline constexpr NodeID InvalidNodeID = AMGraph::InvalidNodeID;
inline constexpr NodeID InvalidEdgeID = AMGraph::InvalidEdgeID;

namespace detail
{
    template <typename T>
    using OrderedSet = Util::BTreeSet<T>;

    template <typename T>
    using UnorderedSet = Util::HashSet<T>;

    template <typename K, typename V, typename Hash>
    using Map = std::unordered_map<K, V, Hash>;

    template <typename K, typename V>
    using OrderedMap = Util::BTreeMap<K, V>;
}


double cal_radians_3d(const Vec3& branch, const Vec3& normal);
double cal_radians_3d(const Vec3& branch, const Vec3& normal, const Vec3& ref_vec);

/// TODO: documentation
struct RsROpts {
    int32_t genus = -1;
    int32_t k = 70;
    double r = 20;
    double theta = 60;
    int32_t n = 50;
    Distance dist = Distance::Euclidean;

    bool is_face_normal = true;
    bool is_face_loop = true;
};

/// Graph definition. The RsR graph here is integrated with the rotation system based on AMGraph
struct Vertex {
    using NormalRep = std::int64_t;

    NormalRep normal_rep = InvalidNormalRep;
    Vec3 coords = Vec3(0., 0., 0.);
    Vec3 normal = Vec3(0., 0., 0.);

    static constexpr NormalRep InvalidNormalRep = -1;
    static constexpr NormalRep CollisionRep = -2;
};

struct Neighbor {
    double angle = 0.0;
    uint v = -1;
    uint tree_id = 0;

    Neighbor(const Vertex& u, const Vertex& v, const uint id)
    {
        this->v = id;
        this->angle = cal_radians_3d(v.coords - u.coords, u.normal);
    }

    // friend size_t hash_value(const Neighbor& p)
    // {
    //     return std::hash<uint>()(p.v);
    // }
    bool operator==(const Neighbor& rhs) const
    {
        return v == rhs.v;
    }

    bool operator<(const Neighbor& rhs) const
    {
        return angle < rhs.angle;
    }

    // std::weak_ordering operator<=>(const Neighbor& rhs) const
    // {
    //     if (this->v == rhs.v) return std::weak_ordering::equivalent;
    //     if (this->angle < rhs.angle) return std::weak_ordering::less;
    //     else return std::weak_ordering::greater;
    // }
};
static_assert(std::is_trivially_destructible_v<Neighbor>);

struct Edge {
    NodeID source = InvalidNodeID;
    NodeID target = InvalidNodeID;
    int ref_time = 0;
};

class SimpGraph /*: public Util::SparseGraph<double>*/ {
    AMGraph graph;
    std::vector<double> m_edges; // Edge weight

public:
    using NodeID = decltype(graph)::NodeID;
    static constexpr auto InvalidEdgeID = decltype(graph)::InvalidEdgeID;

    void reserve(size_t vertices, int k)
    {
        m_edges.reserve(vertices * k);
    }

    // connected components
    [[nodiscard]]
    auto inner() const -> const decltype(graph)&
    {
        return graph;
    }

    // generic algorithms
    auto add_node() -> NodeID
    {
        return graph.add_node();
    }

    auto connect_nodes(const NodeID source, const NodeID target, const double weight = 0.)
    {
        const auto id =
            graph.connect_nodes(source, target);
        if (id == m_edges.size())
            m_edges.push_back(weight);
        else
            m_edges[id] = weight;
        return id;
    }

    [[nodiscard]] double get_weight(const NodeID n1, const NodeID n2) const
    {
        return m_edges[graph.find_edge(n1, n2)];
    }

    /// Disconnect nodes. This operation removes the edge from the edge maps of the two formerly connected
    /// vertices, but the number of edges reported by the super class AMGraph is not decremented, so the edge is only
    /// invalidated. Call cleanup to finalize removal.
    void disconnect_nodes(const NodeID n0, const NodeID n1)
    {
        //return graph.remove_edge(n0, n1);
        //remove_edge(n0, n1);
        // if (graph.valid_node_id(n0) && graph.valid_node_id(n1)) {
        //     graph.edge_map[n0].erase(n1);
        //     graph.edge_map[n1].erase(n0);
        // }
        graph.erase_edge(n0, n1);
    }

    [[nodiscard]]
    auto find_edge(NodeID from, NodeID to) const
    {
        return graph.find_edge(from, to);
    }

    [[nodiscard]]
    auto node_ids() const
    {
        return graph.node_ids();
    }

    [[nodiscard]]
    auto no_nodes() const -> size_t
    {
        return graph.no_nodes();
    }

    [[nodiscard]]
    auto no_edges() const -> size_t
    {
        return graph.no_edges();
    }

    [[nodiscard]]
    auto neighbors_lazy(NodeID from) const
    {
        return graph.neighbors_lazy(from);
    }

    auto collapse(NodeID to_keep, NodeID to_remove)
    {
        for (auto n: neighbors_lazy(to_remove)) {
            double weight = get_weight(n, to_remove);
            connect_nodes(to_keep, n, weight);
        }
        graph.erase_node(to_remove);

    }
};

class RSGraph : private AMGraph {
public:
    static constexpr auto InvalidEdgeID = std::nullopt;

    Geometry::ETF etf;
    std::vector<Vertex> m_vertices;
    std::vector<detail::OrderedSet<Neighbor>> m_neighbors;
    std::vector<Edge> m_edges;

    void reserve(size_t nodes, int k)
    {
        m_vertices.reserve(nodes * k);
    }

    [[nodiscard]]
    auto no_nodes() const -> size_t
    {
        return AMGraph::no_nodes();
    }

    [[nodiscard]]
    auto neighbors_lazy(NodeID n) const
    {
        return AMGraph::neighbors_lazy(n);
    }

    [[nodiscard]]
    auto shared_neighbors_lazy(NodeID n1, NodeID n2) const
    {
        return AMGraph::shared_neighbors_lazy(n1, n2);
    }

    [[nodiscard]]
    auto valence(NodeID n) const
    {
        return AMGraph::valence(n);
    }

    EdgeID add_edge(const NodeID source, const NodeID target)
    {
        const EdgeID id = this->connect_nodes(source, target);
        GEL_ASSERT_NEQ(id, AMGraph::InvalidEdgeID);
        if (m_edges.size() == id)
            m_edges.emplace_back(source, target);
        else
            m_edges[id] = Edge {.source = source, .target = target};

        // insert neighbors
        m_neighbors[source].emplace(Neighbor(m_vertices[source], m_vertices[target], target));
        m_neighbors[target].emplace(Neighbor(m_vertices[target], m_vertices[source], source));

        //return { source, target };
        return id;
    }

    NodeID add_node(const Vec3& p, const Vec3& in_normal = Vec3(0., 0., 0.))
    {
        const NodeID n = AMGraph::add_node();
        GEL_ASSERT_EQ(m_vertices.size(), n);
        m_vertices.emplace_back(Vertex::InvalidNormalRep, p, in_normal);
        m_neighbors.emplace_back();
        //m_vertices[n] = Vertex{.id = n, .normal_rep = Vertex::InvalidNormalRep, .coords = p, .normal = in_normal };
        return n;
    }

    void increment_ref_time(NodeID n1, NodeID n2)
    {
        auto edge = AMGraph::find_edge(n1, n2);
        if (edge != AMGraph::InvalidEdgeID) {
            m_edges[edge].ref_time += 1;
        }
    }

    [[nodiscard]]
    std::optional<Edge> find_edge(NodeID n1, NodeID n2) const
    {
        const EdgeID id = AMGraph::find_edge(n1, n2);
        if (id == AMGraph::InvalidEdgeID) {
            return std::nullopt;
        }
        return m_edges[id];
    }

    /// @brief Get last neighbor information
    /// @param root: root vertex index
    /// @param branch: current outgoing branch
    /// @return reference to last neighbor struct
    // TODO: const correctness
    [[nodiscard]]
    Neighbor& predecessor(const NodeID& root, const NodeID& branch) const
    {
        GEL_ASSERT(m_vertices.size() > root);
        GEL_ASSERT(m_vertices.size() > branch);
        auto& u = m_vertices.at(root);
        auto& v = m_vertices.at(branch);
        GEL_ASSERT(!m_neighbors[root].empty()); // we need at least one neighbor to return
        Neighbor temp = {u, v, static_cast<uint>(branch)};
        auto iter = m_neighbors[root].lower_bound(temp);
        if (iter == m_neighbors[root].begin()) iter = m_neighbors[root].end(); // Wrap around
        //auto& f = iter->first;
        //auto& s = iter->second;
        //return {f, s};
        //return std::tie(f, s);
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
        auto& u = m_vertices.at(root);
        auto& v = m_vertices.at(branch);
        GEL_ASSERT(!m_neighbors[root].empty()); // we need at least one neighbor to return
        auto iter = m_neighbors[root].upper_bound(Neighbor(u, v, branch));
        if (iter == m_neighbors[root].end()) iter = m_neighbors[root].begin(); // Wrap around
        //auto& f = iter->first;
        //auto& s = iter->second;
        //return {f, s};
        //return std::tie(f, s);
        //return const_cast<Neighbor&>(*(iter));
        return (*iter);
    }

private:
    /// @brief Get the neighbor information
    /// @param root: root vertex index
    /// @param branch: the outgoing branch
    /// @return reference to the neighbor struct
    Neighbor&  get_neighbor_info(const NodeID root, const NodeID branch)
    {
        auto& u = this->m_vertices.at(root);
        auto& v = this->m_vertices.at(branch);
        auto iter = m_neighbors[root].lower_bound({u, v, static_cast<uint>(branch)});
        // TODO: tree_id does not invalidate ordered_neighbors, but it still has thread safety issues
        //auto& f = iter->first;
        //auto& s = iter->second;
        //return {f, s};
        return const_cast<Neighbor&>(*(iter));
    }
public:

    void maintain_face_loop(const NodeID source, const NodeID target)
    {
        auto& this_v_tree   = this->predecessor(source, target).tree_id;
        auto& neighbor_tree = this->predecessor(target, source).tree_id;

        auto [fst, snd] = this->etf.insert(this_v_tree, neighbor_tree);

        get_neighbor_info(source, target).tree_id = fst;
        get_neighbor_info(target, source).tree_id = snd;
    }
};

// Algorithm

/// Convert point cloud to a Manifold
/// @param vertices_in vertices of the point cloud
/// @param normals_in normals of the point cloud or empty vector
/// @param opts collapse options
/// @return reconstructed manifold mesh
auto point_cloud_to_mesh(const std::vector<Point>& vertices_in,
                         const std::vector<Vec3>& normals_in,
                         const RsROpts& opts) -> ::HMesh::Manifold;

struct NormalEstimationResult {
    std::vector<Point> vertices;
    std::vector<Vec3> normals;
    std::vector<Vec3> smoothed_v;
};

auto point_cloud_normal_estimate(const std::vector<Point>& vertices,
                                 const std::vector<Vec3>& normals,
                                 bool is_euclidean) -> NormalEstimationResult;


auto point_cloud_collapse_reexpand(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const CollapseOpts& collapse_options,
    const RsROpts& reconstruction_options,
    const ReexpandOptions& reexpand_options) -> HMesh::Manifold;
} // namespace HMesh::RSR

#endif // GEL_HMesh_RsR2_hpp
