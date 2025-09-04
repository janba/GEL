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
#include <unordered_map>
#include <vector>
#include <numbers>

#include <GEL/Geometry/QEM.h>

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

template <std::ranges::input_range R1, std::ranges::input_range R2>
struct zip_view_t : std::ranges::view_interface<zip_view_t<R1, R2>> {
    using value_type = std::pair<std::ranges::range_value_t<R1>, std::ranges::range_value_t<R2>>;
    using reference = value_type;
private:
    std::ranges::iterator_t<R1> m_begin1;
    std::ranges::iterator_t<R2> m_begin2;
    std::ranges::sentinel_t<R1> m_end1;
    std::ranges::sentinel_t<R2> m_end2;
    class iterator_t {
        std::ranges::iterator_t<R1> ptr1;
        std::ranges::iterator_t<R2> ptr2;
    public:
        using difference_type = ptrdiff_t;
        iterator_t() = default;
        explicit iterator_t(auto _ptr1, auto _ptr2) : ptr1(_ptr1), ptr2(_ptr2) {}
        iterator_t& operator++()
        {
            ++ptr1; ++ptr2;
            return *this;
        }
        iterator_t operator++(int) { iterator_t retval = *this; ++(*this); return retval; }
        bool operator==(iterator_t other) const
        {
            return ptr1 == other.ptr1 || ptr2 == other.ptr2;
        }
        bool operator!=(iterator_t other) const { return !(*this == other); }
        reference operator*() const
        {
            return std::make_pair(*ptr1, *ptr2);
        }
    };

public:
    zip_view_t() = default;
    zip_view_t(R1&& range1, R2&& range2) :
        m_begin1(std::ranges::begin(range1)),
        m_begin2(std::ranges::begin(range2)),
        m_end1(std::ranges::end(range1)),
        m_end2(std::ranges::end(range2)) {}

    [[nodiscard]]
    auto begin() const
    {
        return iterator_t(m_begin1, m_begin2);
    }
    [[nodiscard]]
    auto end() const
    {
        return iterator_t(m_end1, m_end2);
    }
};

template <std::ranges::input_range R1, std::ranges::input_range R2>
auto zip(R1&& range1, R2&& range2) -> zip_view_t<R1, R2>
{
    return zip_view_t<R1, R2>{std::forward<R1>(range1), std::forward<R2>(range2)};
}

template <std::ranges::input_range R1, std::ranges::input_range R2>
std::ranges::input_range auto cartesian_product(R1&& range1, R2&& range2)
{
    using T1 = std::ranges::range_value_t<R1>;
    using T2 = std::ranges::range_value_t<R2>;

    auto rv1 = range1 | std::views::all;
    auto rv2 = range2 | std::views::all;

    auto size1 = std::ranges::size(rv1);
    auto size2 = std::ranges::size(rv2);

    auto wide1 = rv1 | std::views::transform([size2](auto e1) {
        std::views::iota(0, size2) | std::views::transform([e1](auto _i) {
            return e1;
        });
    }) | std::views::join;

    auto wide2 = std::views::iota(0, size1) | std::views::transform([rv2](auto _i) {
        return rv2;
    }) | std::views::join;

    return zip(wide1, wide2);
}

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

struct QuadraticCollapseGraph : AMGraph {
    friend struct Collapse;
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

        const auto projected = project_point_to_line(position, m_vertices[latent].position, m_vertices[active].position);
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
        for (const auto n: neighbors_lazy(active)) {
            connect_nodes(n, active);
        }
        edge_map[latent].clear();
        return active;
    }

    struct SingleCollapse {
        NodeID active = InvalidNodeID;
        NodeID latent = InvalidNodeID;
        Point active_point_coords;
        Point latent_point_coords;
        Point v_bar;
    };
    auto collapse_one() -> SingleCollapse
    {
        const auto [edge, _] = m_collapse_queue.pop_front();

        const auto active = m_edges[edge].n0;
        const auto latent = m_edges[edge].n1;
        const auto active_coordinate = m_vertices[active].position;
        const auto latent_coordinate = m_vertices[latent].position;

        const auto v_bar = quadratic_distance(active, latent).first;
        merge_nodes(latent, active, v_bar);

        return { active, latent, active_coordinate, latent_coordinate, v_bar };
    }

    /// Return the remaining points
    auto to_point_cloud() -> PointCloud
    {
        std::vector<Point> points;
        std::vector<Vec3> normals;
        for (auto i = 0UL; i < m_vertices.size(); ++i) {
            if (!m_vertices[i].position.any([](const double e) { return std::isnan(e);})) {
                points.emplace_back(m_vertices[i].position);
                normals.emplace_back(m_vertices[i].normal);
            }
        }
        return PointCloud{std::move(points), std::move(normals)};
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
            const auto v0 = m_vertices[n0];
            const auto v1 = m_vertices[n1];
            const auto qem = v0.qem + v1.qem;
            // We want to clamp the position to be in a sphere bounded by the two points
            // This will get rid of potential outliers
            const auto center = (v0.position + v1.position) * 0.5;
            const auto radius = CGLA::length(v0.position - v1.position) * 0.5;

            const auto shared_neighbors = std::ranges::distance(shared_neighbors_lazy(n0, n1));

            const auto opt_position = qem.opt_pos(0.5, center);
            const auto opt_direction = opt_position - center;
            const auto clamped = center + CGLA::normalize(opt_direction) * std::clamp(CGLA::length(opt_direction), 0.0, radius);
            const auto error = qem.error(clamped);
            // TODO: We can also look for surrounding triangles
            return std::make_pair(clamped, error * radius * radius / (1));
        } else {
            return std::make_pair(CGLA::Vec3d(), CGLA::CGLA_NAN);
        }
    }
};

// Fat 64 bytes
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

// TODO: nuke this class
struct Collapse {

private:
    struct CollapseInfo {
        size_t begin;
        size_t end;
    };

    Util::HashSet<NodeID> m_remaining;
    std::vector<SingleCollapse> m_collapses;
    std::vector<CollapseInfo> m_collapse_ranges;
    PointCloud p;

    Collapse() = default;
public:
    explicit Collapse(std::vector<NodeID>&& remaining): m_remaining(remaining.begin(), remaining.end()) {}

    using ActiveNodeID = NodeID;
    using LatentNodeID = NodeID;
    using CollapseSpan = std::span<const SingleCollapse>;

    struct ActivityMap {
        friend struct Collapse;

    private:
        struct MapVal {
            LatentNodeID latent = InvalidNodeID;
            Point active_point_coords;
            Point latent_point_coords;
            Point v_bar;
        };
        std::vector<std::pair<ActiveNodeID, MapVal>> activity;

    public:
        auto insert(ActiveNodeID active, LatentNodeID latent, const Point& active_point_coords,
                    const Point& latent_point_coords, const Point& v_bar) -> void
        {
            GEL_ASSERT_NEQ(active, latent);
            GEL_ASSERT_NEQ(active, InvalidNodeID);
            GEL_ASSERT_NEQ(latent, InvalidNodeID);
            auto a = std::make_pair(active, MapVal(latent, active_point_coords, latent_point_coords, v_bar));
            activity.emplace_back(a);
        }
    };

    static constexpr NodeID InvalidNodeID = std::numeric_limits<NodeID>::max();

    struct Iterator {
        using difference_type = std::ptrdiff_t;
        using value_type = CollapseSpan;
        friend struct Collapse;
        std::optional<std::reference_wrapper<const Collapse>> collapse;
        size_t element;

        Iterator& operator++()
        {
            element++;
            return *this;
        }

        Iterator operator++(int)
        {
            ++*this;
            return *this;
        }

        Iterator& operator--()
        {
            element--;
            return *this;
        }

        Iterator operator--(int)
        {
            --*this;
            return *this;
        }

        bool operator==(const Iterator& other) const { return element == other.element; }
        bool operator!=(const Iterator& other) const { return element != other.element; }
        CollapseSpan operator*() const { return collapse->get().get_collapse_span(element); }
    };

    static_assert(std::bidirectional_iterator<Iterator>);

    [[nodiscard]] auto begin() const -> Iterator { return Iterator{*this, 0}; }

    [[nodiscard]] auto end() const -> Iterator { return Iterator{*this, number_of_collapses()}; }

    [[nodiscard]]
    auto number_of_collapses() const -> size_t
    {
        return m_collapse_ranges.size();
    }

    // auto insert_one(SingleCollapse collapse) -> void
    // {
    //     m_collapses.emplace_back(collapse);
    // }

    auto insert_collapse(const ActivityMap& activity_map) -> void
    {
        const auto begin_idx = m_collapses.size();
        size_t items = 0;
        for (const auto& info : activity_map.activity | std::views::values) {
            const auto latent = info.latent;
            m_remaining.erase(latent);
            if (latent != InvalidNodeID) {
                const auto active_point_coords = info.active_point_coords;
                const auto latent_point_coords = info.latent_point_coords;
                m_collapses.emplace_back(active_point_coords, latent_point_coords, info.v_bar);
                ++items;
            }
        }
        m_collapse_ranges.emplace_back(CollapseInfo{begin_idx, begin_idx + items});
    }

    [[nodiscard]]
    auto get_collapse_span(const size_t at) const -> CollapseSpan
    {
        if (at >= m_collapse_ranges.size()) { throw std::out_of_range("collapse_ranges"); }
        return {
            m_collapses.begin() + m_collapse_ranges.at(at).begin,
            m_collapses.begin() + m_collapse_ranges.at(at).end
        };
    }

    auto finalize(QuadraticCollapseGraph&& graph) -> void
    {
        auto [points, normals] = std::forward<QuadraticCollapseGraph>(graph).to_point_cloud();
        std::cout << "remaining vertices: " << m_remaining.size() << std::endl;
        p.points = std::move(points);
        p.normals = std::move(normals);
    }

    [[nodiscard]]
    auto remaining() const -> PointCloud
    {
        return p;
    }
};

/// Contains data needed for a reexpansion
struct Collapse2 {
    std::vector<std::vector<SingleCollapse>> collapses;
    // FIXME Should store separately
    //PointCloud remaining;
};

inline auto create_collapse(const Collapse& collapse) -> Collapse2
{
    Collapse2 collapse2;
    for (auto&& iter: collapse) {
        std::vector<SingleCollapse> single_collapses;
        for (auto&& item: iter) {
            single_collapses.push_back(item);
        }
        collapse2.collapses.push_back(std::move(single_collapses));
    }
    return collapse2;
}

inline auto create_collapse(std::vector<SingleCollapse>&& collapses) -> Collapse2
{
    return Collapse2 { .collapses = { std::move(collapses) } };
}

static_assert(std::ranges::viewable_range<Collapse&>);
static_assert(std::ranges::viewable_range<const Collapse&>);

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
            neighbors.reserve(k);
        }
    } else if (neighbors_memoized.at(0).capacity() < k) {
        // TODO: this seems to lose performance because calls to reserve block yet the allocator is multithreaded
        // for (auto &neighbors : neighbors_memoized) {
        //     neighbors.reserve(k);
        // }
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

// TODO: move to cpp file
inline auto collapse_points(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const CollapseOpts& options
) -> Collapse
{
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    Util::ImmediatePool pool;
    QuadraticCollapseGraph graph;

    // initialize graph
    for (auto i = 0UL; i < vertices.size(); ++i) {
        graph.add_node(vertices[i], normals[i]);
    }

    // set of connectivity information
    auto indices = [&vertices] {
        std::vector<NodeID> temp(vertices.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    const auto kd_tree = build_kd_tree_of_indices(vertices, indices);
    const auto neighbor_map = calculate_neighbors(pool, vertices, kd_tree, options.initial_neighbors);
    Collapse collapse{std::move(indices)};

    // This also initializes distances
    for (const auto& neighbors : neighbor_map) {
        const NodeID this_id = neighbors[0].id;
        for (const auto& neighbor : neighbors | std::views::drop(1)) {
            // kNN connection
            graph.connect_nodes(this_id, neighbor.id);
        }
    }

    size_t count = 0;
    for (size_t iter = 0; iter < options.max_iterations; ++iter) {
        // TODO: stricter checking
        const size_t max_collapses = vertices.size() * std::pow(0.5, iter) * options.reduction_per_iteration;
        Collapse::ActivityMap activity_map;

        while (count < max_collapses) {
            count++;

            auto [active, latent, active_point_coords, latent_point_coords, v_bar] = graph.collapse_one();
            
            activity_map.insert(active, latent, active_point_coords, latent_point_coords, v_bar);
        }
        collapse.insert_collapse(activity_map);
    }
    std::cout << "collapsed: " << count << std::endl;
    collapse.finalize(std::move(graph));
    return collapse;
}

struct PointHash {
    size_t operator()(const Point& point) const
    {
        const auto h1 = std::hash<double>{}(point[0]);
        const auto h2 = std::hash<double>{}(point[1]);
        const auto h3 = std::hash<double>{}(point[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct PointEquals {
    size_t operator()(const Point& left, const Point& right) const
    {
        return (left[0] == right[0]) && (left[1] == right[1]) && (left[2] == right[2]);
    }
};

inline Vec3 half_edge_direction(const HMesh::Manifold& m, HMesh::HalfEdgeID h)
{
    const auto w = m.walker(h);
    const auto current = w.vertex();
    const auto opposing = w.opp().vertex();
    return CGLA::normalize(m.positions[opposing] - m.positions[current]);
}

inline Vec3 triangle_normal(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    const auto v1 = p2 - p1;
    const auto v2 = p3 - p1;
    return CGLA::normalize(CGLA::cross(v1, v2));
}

// returns 0 at 180 degrees, 1 at 90 (or 270) degrees
inline double optimize_dihedral(const Vec3& n1, const Vec3& n2)
{
    const auto angle = CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0;
    return std::abs(angle);
    //const auto angle_cos = std::abs(CGLA::dot(n1, n2)) / (CGLA::length(n1) * CGLA::length(n2));
    //GEL_ASSERT_FALSE(std::isnan(angle_cos));
    //const auto angle = std::acos(angle_cos);
    //GEL_ASSERT_FALSE(std::isnan(angle));
    //return angle;
    //return std::abs(angle - 1.0);
}

// returns 0 for an equilateral triangle,
inline double optimize_angle(const Vec3& p1, const Vec3& p2, const Vec3& p3, const ReexpandOptions& opts)
{
    const auto e1 = CGLA::normalize(p2 - p1);
    const auto e2 = CGLA::normalize(p3 - p1);

    const auto e3 = -e2;
    const auto e4 = CGLA::normalize(p3 - p2);

    const auto e5 = -e4;
    const auto e6 = -e1;

    const auto angle1 = std::acos(dot(e1, e2));
    const auto angle2 = std::acos(dot(e3, e4));
    const auto angle3 = std::acos(dot(e5, e6));
    if (angle1 > opts.angle_threshold || angle2 > opts.angle_threshold || angle3 > opts.angle_threshold) {
        return opts.angle_threshold_penalty;
    }

    const auto score1 = std::abs(angle1 - std::numbers::pi / 3.0);
    const auto score2 = std::abs(angle2 - std::numbers::pi / 3.0);
    const auto score3 = std::abs(angle3 - std::numbers::pi / 3.0);
    const auto score = score1 + score2 + score3;
    return score * opts.angle_factor;
}

// penalizes based on effect to valency
// returns 0 if the valency of all affected vertices is 0 as a result of the split
inline double optimize_valency(const HMesh::Manifold& m, const HMesh::HalfEdgeID h_out,
                               const HMesh::HalfEdgeID h_in_opp, const ReexpandOptions& opts)
{
    const auto original_valency = m.valency(m.walker(h_out).opp().vertex());
    // we count the number of edges steps to go from h_in_opp to h_out
    auto steps = 0;
    for (auto w = m.walker(h_out); w.halfedge() != h_in_opp; w = w.circulate_vertex_ccw()) {
        ++steps;
    }
    const auto split_vert_valency = steps + 2;
    const auto original_final_latency = original_valency - steps + 2;

    const auto h_out_end_valency = m.valency(m.walker(h_out).vertex()) + 1;
    const auto h_in_opp_end_valency = m.valency(m.walker(h_in_opp).vertex()) + 1;
    if (split_vert_valency > opts.valency_max_threshold || split_vert_valency < opts.valency_min_threshold
        || original_final_latency > opts.valency_max_threshold || original_final_latency < opts.valency_min_threshold
        || h_out_end_valency > opts.valency_max_threshold || h_out_end_valency < opts.valency_min_threshold
        || h_in_opp_end_valency > opts.valency_max_threshold || h_in_opp_end_valency < opts.valency_min_threshold) {
        return opts.valency_threshold_penalty;
    }

    const auto valency_score =
        std::abs(split_vert_valency - 6)
        + std::abs(original_final_latency - 6)
        + std::abs(h_out_end_valency - 6)
        + std::abs(h_in_opp_end_valency - 6);

    return (valency_score * opts.valency_factor);
}


struct Split {
    HMesh::HalfEdgeID h_in;
    HMesh::HalfEdgeID h_out;
};

struct Triangle {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;

    Vec3 normal() const
    {
        return triangle_normal(p1, p2, p3);
    }
    double area() const
    {
        return triangle_area(p1, p2, p3);
    }
};

inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
    auto e1 = t.p2 - t.p1;
    auto e2 = t.p3 - t.p2;
    auto e3 = t.p1 - t.p3;
    os
        << "{\n"
        << "  p1: " << t.p1 << "\n"
        << "  p2: " << t.p2 << "\n"
        << "  p3: " << t.p3 << "\n"
        << "  e1: " << e1 << "\n"
        << "  e2: " << e2 << "\n"
        << "  e3: " << e3 << "\n"
        << "  normal: " << t.normal() << "\n"
        << "  area: " << t.area() << "\n"
        << "}";
    return os;
}

inline Split find_edge_pair(
    const HMesh::Manifold& m,
    const HMesh::VertexID center_idx,
    const Vec3& v_new_position,
    const Vec3& v_old_position,
    const ReexpandOptions& opts)
{
    // this is getting too complicated

    const auto print_hedge = [&m](HalfEdgeID he) {
        auto v_from = m.positions[m.walker(he).opp().vertex()];
        auto v_to = m.positions[m.walker(he).vertex()];
        std::cout << v_from << " -> " << v_to << " (" << (v_to - v_from) << ")\n";
    };

    const auto triangle_from_half_edge_orig = [&m](const Point& origin, const Walker& w) -> Triangle {
        auto v1 = w.vertex();
        auto v2 = w.next().vertex();
        return Triangle {origin, m.positions[v1], m.positions[v2]};
    };

    const auto triangle_from_half_edge = [&m](const Walker& w) -> Triangle {
        auto v1 = w.vertex();
        auto v2 = w.next().vertex();
        auto v3 = w.next().next().vertex();
        return Triangle {m.positions[v1], m.positions[v2], m.positions[v3]};
    };

    const auto dihedral_ = [](const Triangle& t1, const Triangle& t2) -> double
    {
        const auto d = optimize_dihedral(t1.normal(), t2.normal());
        return d / (t1.area() + t2.area() + DBL_EPSILON);
    };

    const auto dihedral_two_ring = [&m](const Triangle& t, const FaceID f) -> double {
        if (f == InvalidFaceID) return 0.0;
        const auto d = optimize_dihedral(t.normal(), m.normal(f));
        return d / (m.area(f) + t.area() + DBL_EPSILON);
    };

    const auto expand_score = [&](HalfEdgeID h_in_opp, HalfEdgeID h_out, const Point& v_new_position, const Point& v_old_position) -> double {
        const auto walker_out = m.walker(h_out);
        const auto walker_in_opp = m.walker(h_in_opp);

        const auto v_h_out = m.positions[walker_out.vertex()];
        const auto v_h_out_top = m.positions[walker_out.next().vertex()];
        const auto v_h_out_bot = m.positions[walker_out.opp().next().vertex()];

        const auto v_h_in = m.positions[walker_in_opp.vertex()];
        const auto v_h_in_top = m.positions[walker_in_opp.opp().next().vertex()];
        const auto v_h_in_bot = m.positions[walker_in_opp.next().vertex()];

        // Have to check up to 14 triangles

        // instead of doing this unscalable stupidity, let's try to be smart and perform a rotation through all of the
        // affected triangles. We basically perform a sliding window and construct the right triangle by
        // either passing v_new or v_old as the third point. The order of the window also matters a lot

        // center (just 2)
        const auto tri_center_in = Triangle (v_old_position, v_new_position, v_h_in);
        const auto tri_center_out = Triangle (v_old_position, v_h_out, v_new_position);

        // h_in_opp and h_out are unique, making this sound
        std::vector<Triangle> triangles;
        double two_ring = 0.0;

        auto walker_prev = walker_out;
        auto walker_next = walker_prev.circulate_vertex_ccw();
        // from out towards in counterclockwise
        triangles.push_back(tri_center_out);
        while (walker_prev.halfedge() != walker_in_opp.halfedge()) {
            auto p2 = m.positions[walker_prev.vertex()];
            auto p3 = m.positions[walker_next.vertex()];

            // to consider the two ring dihedrals, we need to get the triangles from the opposite edges.
            // none of the triangles are affected by the expansion, so we can just fetch them from the manifold directly
            FaceID opposing_face = walker_next.next().opp().face();
            Triangle tri = {v_new_position, p2, p3};

            two_ring += dihedral_two_ring(tri, opposing_face);

            triangles.push_back(tri);
            walker_prev = walker_prev.circulate_vertex_ccw();
            walker_next = walker_next.circulate_vertex_ccw();
        }
        triangles.push_back(tri_center_in);
        // from in towards out counterclockwise
        walker_prev = walker_in_opp;
        walker_next = walker_in_opp.circulate_vertex_ccw();
        while (walker_prev.halfedge() != walker_out.halfedge()) {
            auto p2 = m.positions[walker_prev.vertex()];
            auto p3 = m.positions[walker_next.vertex()];

            FaceID opposing_face = walker_next.next().opp().face();
            Triangle tri = {v_old_position, p2, p3};

            two_ring += dihedral_two_ring(tri, opposing_face);

            triangles.push_back(tri);
            walker_prev = walker_prev.circulate_vertex_ccw();
            walker_next = walker_next.circulate_vertex_ccw();
        }

        double total_dihedral = 0.0;
        for (size_t i = 0; i < triangles.size(); ++i) {
            const auto& tri1 = triangles[i];
            const auto& tri2 = triangles[(i + 1) % triangles.size()];

            const auto one_ring = dihedral_(tri1, tri2);

            total_dihedral += one_ring;
        }
        total_dihedral += two_ring;

        // top row
        const auto tri_top_in_l = triangle_from_half_edge(walker_in_opp.opp().next().next().opp());
        const auto tri_top_in = Triangle (v_new_position, v_h_in_top, v_h_in);
        const auto tri_top_in_r = triangle_from_half_edge_orig(v_new_position, walker_in_opp.opp().next().opp().next());

        const auto tri_top_out_l = triangle_from_half_edge_orig(v_new_position, walker_out.next().next().opp());
        const auto tri_top_out = Triangle (v_h_out, v_h_out_top, v_new_position);
        const auto tri_top_out_r = triangle_from_half_edge(walker_out.next().opp());

        // bottom row
        const auto tri_bot_in_l = triangle_from_half_edge(walker_in_opp.next().opp());
        const auto tri_bot_in = Triangle (v_old_position, v_h_in, v_h_in_bot);
        const auto tri_bot_in_r = triangle_from_half_edge_orig(v_old_position, walker_in_opp.next().next().opp());

        const auto tri_bot_out_l = triangle_from_half_edge_orig(v_old_position, walker_out.opp().next().opp().next());
        const auto tri_bot_out = Triangle (v_old_position, v_h_out_bot, v_h_out);
        const auto tri_bot_out_r = triangle_from_half_edge(walker_out.opp().next().next().opp());

        auto dihedral0 = dihedral_(tri_center_in, tri_center_out);
        auto dihedral1 = dihedral_(tri_center_in, tri_top_in);
        auto dihedral2 = dihedral_(tri_center_in, tri_bot_in);
        auto dihedral3 = dihedral_(tri_center_out, tri_top_out);
        auto dihedral4 = dihedral_(tri_center_out, tri_bot_out);

        auto dihedral5 = dihedral_(tri_top_in, tri_top_in_r);
        auto dihedral6 = dihedral_(tri_top_out_l, tri_top_out);
        auto dihedral7 = dihedral_(tri_bot_in, tri_bot_in_r); // this seems to reduce performance ?;
        auto dihedral8 = dihedral_(tri_bot_out_l, tri_bot_out);

        if constexpr (DEBUG_PRINT) {
            std::cout << "h_in_opp:\n";
            print_hedge(h_in_opp);
            std::cout << "h_out:\n";
            print_hedge(h_out);
            std::cout << "\n";

            std::cout << "center in: ";
            std::cout << tri_center_in << "\n";
            std::cout << "center out: ";
            std::cout << tri_center_out << "\n";

            std::cout << "top in: ";
            std::cout << tri_top_in << "\n";
            std::cout << "top in r: ";
            std::cout << tri_top_in_r << "\n";

            std::cout << "top out: ";
            std::cout << tri_top_out << "\n";
            std::cout << "top out l: ";
            std::cout << tri_top_out_l << "\n";

            std::cout << "bot in: ";
            std::cout << tri_bot_in << "\n";
            std::cout << "bot in r: ";
            std::cout << tri_bot_in_r << "\n";

            std::cout << "bot out: ";
            std::cout << tri_bot_out << "\n";
            std::cout << "bot out l:";
            std::cout << tri_bot_out_l << "\n";

            std::cout << "center in  <-> center_out: " << dihedral0 << "\n";
            std::cout << "center in  <-> top in:     " << dihedral1 << "\n";
            std::cout << "center in  <-> bot in:     " << dihedral2 << "\n";
            std::cout << "center out <-> top out:    " << dihedral3 << "\n";
            std::cout << "center out <-> bot out:    " << dihedral4 << "\n";
            std::cout << "top in     <-> top in r:   " << dihedral5 << "\n";
            std::cout << "top out l  <-> top out:    " << dihedral6 << "\n";
            std::cout << "bot in     <-> bot in r:   " << dihedral7 << "\n";
            std::cout << "bot out l  <-> bot out:    " << dihedral8 << "\n";
        }

        auto valency_cost = 0; //optimize_valency(m, h_out, h_in_opp, opts);

        return total_dihedral + dihedral0 + valency_cost;
    };

    // let's do it the dumb way for once
    std::vector<HalfEdgeID> half_edges;
    HMesh::circulate_vertex_ccw(m, center_idx, [&](HalfEdgeID h) {
       half_edges.push_back(h);
    });

    if constexpr (DEBUG_PRINT) {
        std::cout << "half_edges: " << half_edges.size() << std::endl;
        std::cout << "v_old: " << v_old_position << "\n";
        std::cout << "v_new: " << v_new_position << "\n";
    }

    double min_score = INFINITY;
    HalfEdgeID h_in_opp;
    HalfEdgeID h_out;
    for (auto h1: half_edges) {
        for (auto h2: half_edges) {
            if (h1 == h2) {
                continue;
            }
            double score = expand_score(h1, h2, v_new_position, v_old_position);
            if (score < min_score) {
                min_score = score;
                h_in_opp = h1;
                h_out = h2;
                // FIXME DEBUG
                if constexpr (DEBUG_PRINT) {
                    std::cout << "<<>>score: " << score << "\n";
                    std::cout << "<<>>h_in_opp: ";
                    print_hedge(h1);
                    std::cout << "<<>>h_out: ";
                    print_hedge(h2);
                    std::cout << "\n";
                }

            } else {
                if constexpr (DEBUG_PRINT) {
                    std::cout << "score: " << score << "\n";
                    std::cout << "h_in_opp: ";
                    print_hedge(h1);
                    std::cout << "h_out: ";
                    print_hedge(h2);
                    std::cout << "\n";
                }
            }
            if constexpr (DEBUG_PRINT) {
                std::cout << "-------------------------" << "\n";
            }
        }
    }

    return Split{m.walker(h_in_opp).opp().halfedge(), h_out};
}


inline auto reexpand_points(::HMesh::Manifold& manifold, Collapse2&& collapse,
                            const ReexpandOptions& opts = ReexpandOptions()) -> void
{
    std::cout << "reexpanding\n";
    const auto& manifold_positions = manifold.positions;

    // TODO: Replace this with the collection in Util
    std::unordered_multimap<Point, VertexID, PointHash, PointEquals> point_to_manifold_ids;
    for (auto manifold_vid : manifold.vertices()) {
        auto pos = manifold_positions[manifold_vid];
        point_to_manifold_ids.emplace(pos, manifold_vid);
    }
    auto position_to_manifold_iter = [&](const Point& point) {
        auto [fst, snd] = point_to_manifold_ids.equal_range(point);
        return std::ranges::subrange(fst, snd) | std::views::values;
    };
    using IndexIter = decltype(position_to_manifold_iter(std::declval<Point>()));

    // insert latent point to stored latent position
    // update active point position to the stored coordinate

    // Now we need to consider two position candidates
    auto try_two_edge_expand = [&](IndexIter vs, const Point& latent_pos, const Point& active_pos) -> VertexID {
        // we want to get as close to 90 degrees as possible here
        for (const auto this_vert : vs) {
            const auto candidate = find_edge_pair(manifold, this_vert, latent_pos, active_pos, opts);
            if (candidate.h_in != InvalidHalfEdgeID && candidate.h_out != InvalidHalfEdgeID) {
                const auto vnew = manifold.split_vertex(candidate.h_in, candidate.h_out);
                GEL_ASSERT_NEQ(vnew, InvalidVertexID);
                return vnew;
            }
        }
        return InvalidVertexID;
    };

    size_t failures = 0;
    for (const auto& collapse_iter : collapse.collapses | std::views::reverse) {
        for (auto single_collapse : collapse_iter | std::views::reverse) {
            // find the manifold_ids for the active vertex
            //const auto active_idx = single_collapse.active;
            //const auto latent_idx = single_collapse.latent;
            const auto active_pos = single_collapse.active_point_coords;
            const auto latent_pos = single_collapse.latent_point_coords;
            const auto v_bar = single_collapse.v_bar;
            // need old active coordinates
            const auto manifold_ids = position_to_manifold_iter(v_bar);
            for (const auto id: manifold_ids) {
                manifold.positions[id] = active_pos;
            }
            if (const auto new_vid = try_two_edge_expand(manifold_ids, latent_pos, active_pos); new_vid != InvalidVertexID) {
                manifold.positions[new_vid] = latent_pos;
                for (const auto id: manifold_ids) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
                point_to_manifold_ids.emplace(latent_pos, new_vid);
                point_to_manifold_ids.erase(v_bar);
            } else {
                failures++;
            }
        }
    }
    std::cerr << "failures: " << failures << "\n";
}

inline auto decimate(const Manifold& manifold, double factor = 0.1) -> Manifold
{
    if (factor >= 1.0 || factor < 0.0) {
        throw std::runtime_error("Invalid factor");
    }

    QuadraticCollapseGraph graph;
    // insert QEM for every point
    for (const auto vertex_id: manifold.vertices()) {
        graph.add_node(manifold.positions[vertex_id], manifold.normal(vertex_id));
    }
    // insert every edge into a queue
    for (const auto edge_id: manifold.halfedges()) {
        auto walker = manifold.walker(edge_id);
        auto p1 = walker.vertex();
        auto p2 = walker.opp().vertex();
        if (p1.get_index() > p2.get_index()) {
            graph.connect_nodes(p1.get_index(), p2.get_index());
        }
    }

    // perform a collapse until we reach the desired number of points
    const size_t max_collapses = manifold.no_vertices() * (1.0 - factor);
    for (size_t i = 0; i < max_collapses; ++i) {
        [[maybe_unused]]
        auto collapse = graph.collapse_one();
    }

    // create a new manifold from the collapsed graph
    std::vector<Vec3> vertices;
    vertices.reserve(graph.no_nodes());
    std::vector<NodeID> indices;
    indices.reserve(graph.no_edges() * 2);
    for (NodeID id: graph.node_ids()) {
        auto p1 = graph.m_vertices.at(id).position;
        vertices.push_back(p1);
        for (NodeID neighbor: graph.neighbors_lazy(id)) {
            if (id < neighbor) {
                for (NodeID third: graph.shared_neighbors(id, neighbor)) {
                    if (neighbor < third) {
                        auto n1 = graph.m_vertices.at(id).normal;
                        auto n2 = graph.m_vertices.at(neighbor).normal;
                        auto n3 = graph.m_vertices.at(third).normal;

                        auto p2 = graph.m_vertices.at(neighbor).position;
                        auto p3 = graph.m_vertices.at(third).position;
                        auto e1 = p2 - p1;
                        auto e2 = p3 - p1;
                        auto n = CGLA::cross(e1, e2);
                        if (CGLA::dot((n1 + n2 + n3) * (1./3.), n) > 0.0) {
                            indices.push_back(id);
                            indices.push_back(neighbor);
                            indices.push_back(third);
                        } else {
                            indices.push_back(id);
                            indices.push_back(third);
                            indices.push_back(neighbor);
                        }
                    }
                }
            }
        }
    }
    Manifold m;
    build_manifold(m, vertices, indices, 3);
    return m;
}

inline auto decimate_reexpand(const Manifold& manifold, double factor = 0.1) -> Manifold
{
        if (factor >= 1.0 || factor < 0.0) {
        throw std::runtime_error("Invalid factor");
    }

    QuadraticCollapseGraph graph;
    // insert QEM for every point
    for (const auto vertex_id: manifold.vertices()) {
        graph.add_node(manifold.positions[vertex_id], manifold.normal(vertex_id));
    }
    // insert every edge into a queue
    for (const auto edge_id: manifold.halfedges()) {
        auto walker = manifold.walker(edge_id);
        auto p1 = walker.vertex();
        auto p2 = walker.opp().vertex();
        if (p1.get_index() > p2.get_index()) {
            graph.connect_nodes(p1.get_index(), p2.get_index());
        }
    }

    Collapse2 collapse;
    collapse.collapses.emplace_back();
    // perform a collapse until we reach the desired number of points
    const size_t max_collapses = manifold.no_vertices() * (1.0 - factor);
    for (size_t i = 0; i < max_collapses; ++i) {
        auto single_collapse = graph.collapse_one();
        collapse.collapses[0].emplace_back(single_collapse.active_point_coords, single_collapse.latent_point_coords, single_collapse.v_bar);
    }

    // create a new manifold from the collapsed graph
    std::vector<Vec3> vertices;
    vertices.reserve(graph.no_nodes());
    std::vector<NodeID> indices;
    indices.reserve(graph.no_edges() * 2);
    for (NodeID id: graph.node_ids()) {
        auto p1 = graph.m_vertices.at(id).position;
        vertices.push_back(p1);
        for (NodeID neighbor: graph.neighbors_lazy(id)) {
            if (id < neighbor) {
                for (NodeID third: graph.shared_neighbors(id, neighbor)) {
                    if (neighbor < third) {
                        auto n1 = graph.m_vertices.at(id).normal;
                        auto n2 = graph.m_vertices.at(neighbor).normal;
                        auto n3 = graph.m_vertices.at(third).normal;

                        auto p2 = graph.m_vertices.at(neighbor).position;
                        auto p3 = graph.m_vertices.at(third).position;
                        auto e1 = p2 - p1;
                        auto e2 = p3 - p1;
                        auto n = CGLA::cross(e1, e2);
                        if (CGLA::dot((n1 + n2 + n3) * (1./3.), n) > 0.0) {
                            indices.push_back(id);
                            indices.push_back(neighbor);
                            indices.push_back(third);
                        } else {
                            indices.push_back(id);
                            indices.push_back(third);
                            indices.push_back(neighbor);
                        }
                    }
                }
            }
        }
    }
    Manifold m;
    build_manifold(m, vertices, indices, 3);

    reexpand_points(m, std::move(collapse));

    return m;
}


} // namespace HMesh

#endif // GEL_HMESH_COLLAPSE_H
