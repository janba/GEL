//
// Created by Cem Akarsubasi on 7/11/25.
// Some Collapse data structures and functions that can be separated from RsR

#ifndef GEL_HMESH_COLLAPSE_H
#define GEL_HMESH_COLLAPSE_H

#include <GEL/HMesh/Manifold.h>
#include <GEL/Util/Assert.h>
#include <GEL/Util/ParallelAdapters.h>

#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>

#include <span>
#include <vector>
#include <numbers>

#include <GEL/Geometry/QEM.h>

#include <GEL/Util/InplaceVector.h>


namespace HMesh::RSR
{
using Vec3 = CGLA::Vec3d;
using Point = Vec3;
using NodeID = size_t;
using Geometry::AMGraph;
using Geometry::AMGraph3D;

static constexpr bool DEBUG_PRINT = false;

struct PointCloud {
    std::vector<Point> points;
    std::vector<Vec3> normals;
    std::vector<NodeID> indices;
};

struct CollapseOpts {
    size_t max_iterations = 1;
    double reduction_per_iteration = 0.5;
    size_t initial_neighbors = 8;
};

struct ReexpandOptions {
    double angle_factor = 0.008 * 0.0; // TODO
    double valency_factor = 0.25 * 0.0; // TODO
    /// maximum angle to allow a triangle
    double angle_threshold = std::numbers::pi * 0.90; // 162 degrees

    double valency_min_threshold = 3.0;
    // double boundary_valency_min_threshold = 4.0;
    double valency_max_threshold = 8.0;
    // double boundary_valency_max_threshold = 8.0;
    // TODO: might be worth considering whether to treat boundaries differently

    double angle_threshold_penalty = 500.0;
    double valency_threshold_penalty = 1000.0;
};

inline auto project_point_to_line(const Vec3& to_project, const Point& p1, const Point& p2) -> double
{
    const auto e1 = p2 - p1;
    const auto e2 = to_project - p1;
    return (CGLA::dot(e1, e2) / CGLA::length(e1));
    //* e1 + p1;
}

inline auto lerp(const Vec3& v1, const Vec3& v2, double t) -> Vec3
{
    return v1 * (1.0 - t) + v2 * t;
}

inline double triangle_area(const Point& p1, const Point& p2, const Point& p3)
{
    const auto e1 = p2 - p1;
    const auto e2 = p3 - p1;
    return CGLA::length(CGLA::cross(e1, e2)) / 2.0;
}

// Fat 72 bytes
struct SingleCollapse {
    //NodeID active = InvalidNodeID;
    //NodeID latent = InvalidNodeID;
    /// Old coordinates of the active point
    Point active_point_coords;
    /// Old coordinates of the latent point
    Point latent_point_coords;
    /// Current coordinates of the active point
    Point v_bar;
};

struct QuadraticCollapseGraph : AMGraph {

private:
    struct Vertex {
        Point position;
        Vec3 normal;
        GEO::QEM qem;
    };

    struct Edge {
        NodeID n0;
        NodeID n1;
    };

public:
    Util::AttribVec<EdgeID, Edge> m_edges;
    Util::AttribVec<NodeID, Vertex> m_vertices;
    Util::MutablePriorityQueue<EdgeID, double, std::hash<EdgeID>, std::equal_to<>, std::less<>> m_collapse_queue;

public:
    auto add_node(const Point& position, const Vec3& normal) -> NodeID
    {
        const NodeID n = AMGraph::add_node();
        const auto qem = Geometry::QEM(position, normal);
        m_vertices[n] = Vertex{.position = position, .normal = normal, .qem = qem};
        return n;
    }

    auto connect_nodes(const NodeID n0, const NodeID n1) -> EdgeID
    {
        const EdgeID e = AMGraph::connect_nodes(n0, n1);
        m_edges[e] = Edge{.n0 = n0, .n1 = n1};
        auto [_, dist] = quadratic_distance(n0, n1);
        std::array elem = {std::make_pair(e, dist)};
        m_collapse_queue.insert_range(elem);

        return e;
    }

    auto merge_nodes(const NodeID latent, const NodeID active, const Point& position) -> NodeID
    {
        GEL_ASSERT_NEQ(latent, active);
        GEL_ASSERT_NEQ(AMGraph::find_edge(latent, active), InvalidEdgeID);

        const auto projected =
            project_point_to_line(position, m_vertices[latent].position, m_vertices[active].position);
        const auto new_normal = lerp(m_vertices[latent].normal, m_vertices[active].normal, projected);
        m_vertices[active].position = position;
        m_vertices[active].normal = new_normal;
        //m_vertices[active].normal = m_vertices[latent].normal;
        m_vertices[active].qem += m_vertices[latent].qem;
        m_vertices[latent].position = Point(std::numeric_limits<double>::signaling_NaN());
        m_vertices[latent].normal = Vec3(std::numeric_limits<double>::signaling_NaN());

        const auto edges_latent = AMGraph::edges(latent) | std::views::values;
        const auto edges_active = AMGraph::edges(active) | std::views::values;
        // clear up old ranges
        m_collapse_queue.remove_range(edges_latent);
        m_collapse_queue.remove_range(edges_active);

        for (const auto n : neighbors_lazy(latent)) {
            edge_map[n].erase(latent);
            if (n != active && n != latent)
                // reconnections handle the new insertions
                connect_nodes(n, active);
        }
        for (const auto n : neighbors_lazy(active)) {
            connect_nodes(n, active);
        }
        edge_map[latent].clear();
        return active;
    }

    struct RawCollapse {
        NodeID active = InvalidNodeID;
        NodeID latent = InvalidNodeID;
        Point active_point_coords;
        Point latent_point_coords;
        Point v_bar;
    };

    auto collapse_one() -> RawCollapse
    {
        const auto [edge, _] = m_collapse_queue.pop_front();

        const auto active = m_edges[edge].n0;
        const auto latent = m_edges[edge].n1;
        const auto active_coordinate = m_vertices[active].position;
        const auto latent_coordinate = m_vertices[latent].position;

        const auto v_bar = quadratic_distance(active, latent).first;
        merge_nodes(latent, active, v_bar);

        return { active, latent, active_coordinate, latent_coordinate, v_bar};
    }

    /// Return the remaining points
    auto to_point_cloud() -> PointCloud
    {
        std::vector<Point> points;
        std::vector<Vec3> normals;
        std::vector<NodeID> indices;
        for (auto i = 0UL; i < m_vertices.size(); ++i) {
            if (!m_vertices[i].position.any([](const double e) { return std::isnan(e); })) {
                points.emplace_back(m_vertices[i].position);
                normals.emplace_back(m_vertices[i].normal);
                indices.emplace_back(i);
            }
        }
        return PointCloud{std::move(points), std::move(normals), std::move(indices)};
    }

private:
    /// Compute sqr distance between two nodes - not necessarily connected.
    [[nodiscard]]
    double sqr_dist(const NodeID n0, const NodeID n1) const
    {
        if (valid_node_id(n0) && valid_node_id(n1))
            return CGLA::sqr_length(m_vertices[n0].position - m_vertices[n1].position);
        else
            return CGLA::CGLA_NAN;
    }

    [[nodiscard]]
    std::pair<CGLA::Vec3d, double> quadratic_distance(const NodeID n0, const NodeID n1) const
    {
        GEL_ASSERT(valid_node_id(n0) && valid_node_id(n1));
        if (valid_node_id(n0) && valid_node_id(n1)) {
            const auto& v0 = m_vertices[n0];
            const auto& v1 = m_vertices[n1];
            const auto qem = v0.qem + v1.qem;
            // We want to clamp the position to be in a sphere bounded by the two points
            // This will get rid of potential outliers
            const auto center = (v0.position + v1.position) * 0.5;
            const auto radius = CGLA::length(v0.position - v1.position) * 0.5;

            const auto opt_position = qem.opt_pos(0.5, center);
            const auto opt_direction = opt_position - center;
            const auto clamped = center + CGLA::normalize(opt_direction) * std::clamp(
                CGLA::length(opt_direction), 0.0, radius);
            const auto error = qem.error(clamped);
            return std::make_pair(clamped, error * radius * radius);
        } else {
            return std::make_pair(CGLA::Vec3d(), CGLA::CGLA_NAN);
        }
    }
};

/// Contains data needed for a reexpansion
struct Collapse {
    std::vector<std::vector<SingleCollapse>> collapses;
};

inline auto create_collapse(std::vector<SingleCollapse>&& collapses) -> Collapse
{
    return Collapse{.collapses = {std::move(collapses)}};
}

using Tree = Geometry::KDTree<Point, NodeID>;
using Record = Geometry::KDTreeRecord<Point, NodeID>;

struct NeighborInfo {
    NodeID id;
    double distance; // the added precision doesn't justify the usage of doubles

    static constexpr NodeID invalid_id = -1;
    //NeighborInfo() = default;

    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
};

//using NeighborArray = Util::InplaceVector<NeighborInfo, 192>;
using NeighborArray = std::vector<NeighborInfo>;

using NeighborMap = std::vector<NeighborArray>;

template <typename Indices>
Tree build_kd_tree_of_indices(const std::vector<Point>& vertices, const Indices& indices)
{
    Tree kd_tree;
    for (const auto idx : indices) {
        kd_tree.insert(vertices.at(idx), idx);
    }
    kd_tree.build();
    return kd_tree;
}

/// @brief k nearest neighbor search
/// @param query: the coordinate of the point to be queried
/// @param kdTree: kd-tree for knn query
/// @param num: number of nearest neighbors to be queried
/// @param neighbors: [OUT] indices of k nearest neighbors
inline void knn_search(const Point& query, const Tree& kdTree, const int num, NeighborArray& neighbors)
{
    // It might be a better idea for the caller to handle this to reduce some clutter
    std::vector<Record> records = kdTree.m_closest(num + 1, query, INFINITY);
    std::sort_heap(records.begin(), records.end());

    for (auto record : records) {
        neighbors.emplace_back(record);
    }
}

template <typename Indices>
auto calculate_neighbors(
    Util::IExecutor& pool,
    const std::vector<Point>& vertices,
    const Indices& indices,
    const Tree& kdTree,
    const int k,
    NeighborMap&& neighbors_memoized = NeighborMap())
    -> NeighborMap
{
    if (neighbors_memoized.empty()) {
        neighbors_memoized = NeighborMap(indices.size());
        for (auto& neighbors : neighbors_memoized) {
            neighbors.reserve(k + 2);
        }
    } else if (neighbors_memoized.at(0).capacity() < k) {
        for (auto &neighbors : neighbors_memoized) {
            neighbors.reserve(k + 2);
        }
    }

    auto cache_kNN_search = [&kdTree, k, &vertices](auto index, auto& neighbor) {
        auto vertex = vertices.at(index);
        knn_search(vertex, kdTree, k, neighbor);
    };
    Util::Parallel::foreach2(pool, indices, neighbors_memoized, cache_kNN_search);
    return neighbors_memoized;
}

inline auto calculate_neighbors(
    Util::IExecutor& pool,
    const std::vector<Point>& vertices,
    const Tree& kdTree,
    const int k,
    NeighborMap&& neighbors_memoized = NeighborMap())
    -> NeighborMap
{
    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    return calculate_neighbors(pool, vertices, indices, kdTree, k, std::move(neighbors_memoized));
}

auto collapse_points(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const CollapseOpts& options
) -> std::pair<Collapse, PointCloud>;

auto reexpand_points(Manifold& manifold, Collapse&& collapse,
                     const ReexpandOptions& opts = ReexpandOptions()) -> void;

auto decimate(const Manifold& manifold, double factor = 0.1) -> Manifold;

auto decimate_reexpand(const Manifold& manifold, double factor = 0.1) -> Manifold;
} // namespace HMesh

#endif // GEL_HMESH_COLLAPSE_H
