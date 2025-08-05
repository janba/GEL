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

#include "GEL/CGLA/Mat4x4d.h"

namespace HMesh::RSR
{
using Vec3 = CGLA::Vec3d;
using Point = Vec3;
using NodeID = size_t;
using Geometry::AMGraph;
using Geometry::AMGraph3D;

struct ReexpandOptions {
    double angle_factor = 0.008 * 0.0; // TODO
    double valency_factor = 0.25 * 0.0; // TODO
    /// maximum angle to allow a triangle
    double angle_threshold = std::numbers::pi * 0.90; // 162 degrees

    double valency_min_threshold = 4.0;
    // double boundary_valency_min_threshold = 4.0;
    double valency_max_threshold = 8.0;
    // double boundary_valency_max_threshold = 8.0;
    // TODO: might be worth considering whether to treat boundaries differently

    // TODO: We end up dealing with 5 triangles for dihedral angles and it might be worth considering bias

    double angle_threshold_penalty = 500.0;
    double valency_threshold_penalty = 1000.0;
};

struct Collapse {
    std::vector<NodeID> m_remaining;

private:
    std::vector<std::pair<NodeID, NodeID>> m_collapses;

    struct CollapseInfo {
        size_t begin;
        size_t end;
    };

    std::vector<CollapseInfo> m_collapse_ranges;

public:
    explicit Collapse(std::vector<NodeID>&& remaining): m_remaining(std::move(remaining)) {}

    using ActiveNodeID = NodeID;
    using LatentNodeID = NodeID;
    using CollapseSpan = std::span<const std::pair<NodeID, NodeID>>;

    struct ActivityMap {
        friend struct Collapse;

    private:
        std::unordered_map<ActiveNodeID, std::pair<LatentNodeID, Point>> activity;

    public:
        auto insert(ActiveNodeID active, LatentNodeID latent, const Point& latent_point_coords) -> void
        {
            GEL_ASSERT_NEQ(active, latent);
            GEL_ASSERT_NEQ(active, InvalidNodeID);
            GEL_ASSERT_NEQ(latent, InvalidNodeID);
            GEL_ASSERT(!activity.contains(active));
            GEL_ASSERT(!activity.contains(latent));
            activity.emplace(active, std::make_pair(latent, latent_point_coords));
            activity.emplace(latent, std::make_pair(Collapse::InvalidNodeID, Point()));
        }

        auto is_used(const NodeID id) const -> bool { return activity.contains(id); }

        auto is_active(const ActiveNodeID id) const -> bool
        {
            if (const auto maybe = activity.find(id); maybe != activity.end()) {
                return maybe->second.first != Collapse::InvalidNodeID;
            }
            return false;
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

    auto insert_collapse(const ActivityMap& activity_map) -> void
    {
        const auto begin_idx = m_collapses.size();
        // const auto end_idx = begin_idx + activity_map.activity.size();
        size_t items = 0;
        for (auto [active, latent] : activity_map.activity) {
            if (latent.first != InvalidNodeID) {
                m_collapses.emplace_back(active, latent.first);
                ++items;
            }
        }
        m_collapse_ranges.emplace_back(CollapseInfo{begin_idx, begin_idx + items});

        // update remaining indices
        auto [fst, snd] = std::ranges::remove_if(m_remaining, [&](const NodeID id) {
            const auto maybe = activity_map.activity.find(id);
            return (maybe != activity_map.activity.end() && maybe->second.first == Collapse::InvalidNodeID);
        });
        m_remaining.erase(fst, snd);
    }

    [[nodiscard]]
    auto get_collapse_span(const size_t at) const -> std::span<const std::pair<NodeID, NodeID>>
    {
        if (at >= m_collapse_ranges.size()) { throw std::out_of_range("collapse_ranges"); }
        return {
            m_collapses.begin() + m_collapse_ranges.at(at).begin,
            m_collapses.begin() + m_collapse_ranges.at(at).end
        };
    }
};

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

/**
    * @brief k nearest neighbor search
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param num: number of nearest neighbors to be queried
    * @param neighbors: [OUT] indices of k nearest neighbors
    *
    * @return None
    */
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

inline auto quadratic_error(const CGLA::Mat4x4d& qem) -> std::pair<CGLA::Vec3d, double>
{
    auto q = CGLA::Mat4x4d(
        qem[0],
        qem[1],
        qem[2],
        CGLA::Vec4d(0.0, 0.0, 0.0, 1.0));
    auto det = CGLA::determinant(q);
    if (CGLA::is_tiny(det)) {
        GEL_ASSERT(false, "shit");
        return std::make_pair(CGLA::Vec3d(), 0.0);
    } else {
        auto v_bar = CGLA::invert(q) * CGLA::Vec4d(0.0, 0.0, 0.0, 1.0);
        CGLA::Vec4d temp;
        for (int i = 0; i < 4; ++i) {
            auto column = CGLA::Vec4d(qem[i][0], qem[i][1], qem[i][2], qem[i][3]);
            temp[i] = CGLA::dot(column, v_bar);
        }
        auto score = CGLA::dot(temp, v_bar);
        return std::make_pair(CGLA::Vec3d(v_bar), score);
    }
}

struct QuadraticCollapseGraph : AMGraph {
    struct Vertex {
        Point position;
        Vec3 normal;
        CGLA::Mat4x4d qem;
    };

    Util::AttribVec<NodeID, Vertex> m_vertices;

    auto add_node(const Point& position, const Vec3& normal) -> NodeID
    {
        const NodeID n = AMGraph::add_node();
        m_vertices[n] = Vertex{.position = position, .normal = normal};
        return n;
    }

    auto calculate_qem() -> void
    {
        for (auto i = 0UL; i < m_vertices.size(); ++i) {
            auto& vertex = m_vertices[i];
            auto p = CGLA::Vec4d(vertex.normal, CGLA::dot(-vertex.normal, vertex.position));
            CGLA::outer_product(p, p, vertex.qem);
        }
    }

    auto merge_nodes(const NodeID n0, const NodeID n1, const Point& position) -> NodeID
    {
        if (n0 == n1) {
            edge_map[n1].erase(n1);
        } else {
            m_vertices[n0].position = Point(std::numeric_limits<double>::quiet_NaN());
            m_vertices[n1].position = position;
            m_vertices[n1].normal = m_vertices[n0].normal;
            m_vertices[n1].qem += m_vertices[n0].qem;

            for (const auto n : neighbors(n0)) {
                edge_map[n].erase(n0);
                if (n != n1 && n != n0)
                    connect_nodes(n, n1);
            }
            edge_map[n0].clear();
        }
        return n1;
    }

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
    std::pair<CGLA::Vec3d, double> quadratic_distance(NodeID n0, NodeID n1) const
    {
        if (valid_node_id(n0) && valid_node_id(n1)) {
            return quadratic_error(m_vertices[n0].qem + m_vertices[n1].qem);
        } else {
            return std::make_pair(CGLA::Vec3d(), CGLA::CGLA_NAN);
        }
    }
};

// TODO: move to cpp file
inline auto collapse_points(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const size_t max_iterations
) -> Collapse
{
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    Util::ImmediatePool pool;
    QuadraticCollapseGraph graph;

    for (auto i = 0UL; i < vertices.size(); ++i) {
        graph.add_node(vertices[i], normals[i]);
    }

    auto indices = [&vertices] {
        std::vector<NodeID> temp(vertices.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    const auto kd_tree = build_kd_tree_of_indices(vertices, indices);
    const auto neighbor_map = calculate_neighbors(pool, vertices, kd_tree, 9);
    Collapse collapse{std::move(indices)};

    for (auto& neighbors : neighbor_map) {
        const NodeID this_id = neighbors[0].id;
        for (auto& neighbor : std::ranges::subrange(neighbors.begin() + 1, neighbors.end())) {
            // kNN connection
            graph.connect_nodes(this_id, neighbor.id);
        }
    }

    struct EdgeElem {
        NodeID n0;
        NodeID n1;
        Point v_bar;
        double dist;

        bool operator<(const EdgeElem& other) const
        {
            return dist < other.dist;
        }
    };
    static_assert(std::movable<EdgeElem>);
    std::priority_queue<EdgeElem> queue;
    auto priority = [&](NodeID a, NodeID b) { return graph.quadratic_distance(a, b); };

    size_t count = 0;
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        Collapse::ActivityMap activity_map;

        Util::AttribVec<NodeID, uint8_t> visited(graph.no_nodes(), 0);
        for (auto n0 : graph.node_ids()) {
            for (auto n1 : graph.neighbors(n0)) {
                auto [v_bar, pri] = priority(n0, n1);
                // TODO: It seems that Util::Range returns a i64 regardless which is probably not what we want
                // TODO: Since I will most likely remove that in exchange for iota_view (despite the apparent downsides)
                // TODO: This will likely not be a problem later.
                queue.emplace(EdgeElem{.n0 = static_cast<NodeID>(n0), .n1 = n1, .v_bar = v_bar, .dist = -pri});
            }
        }
        while (!queue.empty()) {
            auto shortest_edge = queue.top();
            queue.pop();

            if (visited[shortest_edge.n0] || visited[shortest_edge.n1]) continue;

            visited[shortest_edge.n0] = 1;
            visited[shortest_edge.n1] = 1;
            count++;
            // TODO: need to store both coordinates
            //const auto old_coordinate
            graph.merge_nodes(shortest_edge.n0, shortest_edge.n1, shortest_edge.v_bar);
            activity_map.insert(shortest_edge.n1, shortest_edge.n0, Point());
        }
        collapse.insert_collapse(activity_map);
    }
    std::cout << "collapsed: " << count << std::endl;

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

inline std::array<HMesh::HalfEdgeID, 2> find_edge_pair(const HMesh::Manifold& m, const HMesh::VertexID center_idx,
                                                       const Vec3& to_insert_position, const ReexpandOptions& opts)
{
    std::array edges = {HMesh::InvalidHalfEdgeID, HMesh::InvalidHalfEdgeID};
    const auto center_position = m.positions[center_idx];

    struct Candidate {
        HalfEdgeID h;
        Vec3 normal;
        double score;
        Vec3 points;

        std::weak_ordering operator<=>(const Candidate& other) const
        {
            return this->score == other.score
                       ? std::weak_ordering::equivalent
                       : this->score < other.score
                       ? std::weak_ordering::less
                       : std::weak_ordering::greater;
        }
    };
    struct CandidatePair {
        HalfEdgeID in;
        HalfEdgeID out;
        double score;

        std::weak_ordering operator<=>(const CandidatePair& other) const
        {
            return this->score == other.score
                       ? std::weak_ordering::equivalent
                       : this->score < other.score
                       ? std::weak_ordering::less
                       : std::weak_ordering::greater;
        }
    };

    std::vector<Candidate> h_out_candidates;
    std::vector<Candidate> h_in_candidates;
    std::priority_queue<CandidatePair, std::vector<CandidatePair>, std::greater<>> h_pair_candidates;
    HMesh::circulate_vertex_ccw(m, center_idx, [&](HMesh::Walker& h) {
        const auto v_h_out = m.positions[h.vertex()];
        const auto v_h_out_top = m.positions[h.next().vertex()];
        const auto v_h_out_bot = m.positions[h.opp().next().vertex()];
        const auto v_h_in = m.positions[h.vertex()];
        const auto v_h_in_top = m.positions[h.opp().next().vertex()];
        const auto v_h_in_bot = m.positions[h.next().vertex()];

        const auto tri_out_top = triangle_normal(to_insert_position, v_h_out_top, v_h_out);
        const auto tri_out_new = triangle_normal(center_position, to_insert_position, v_h_out);
        const auto tri_out_bot = triangle_normal(v_h_out, v_h_out_bot, center_position);
        const auto out_dihedral_score
            = optimize_dihedral(tri_out_top, tri_out_new) + optimize_dihedral(tri_out_new, tri_out_bot);

        const auto sharpness_tri_out_top = optimize_angle(to_insert_position, v_h_out_top, v_h_out, opts);
        const auto sharpness_tri_out_new = optimize_angle(center_position, to_insert_position, v_h_out, opts);
        const auto sharpness_tri_out_bot = optimize_angle(v_h_out, v_h_out_bot, center_position, opts);
        const auto out_sharpness = sharpness_tri_out_bot + sharpness_tri_out_new + sharpness_tri_out_top;

        const auto tri_in_top = triangle_normal(to_insert_position, v_h_in, v_h_in_top);
        const auto tri_in_new = triangle_normal(to_insert_position, center_position, v_h_in);
        const auto tri_in_bot = triangle_normal(center_position, v_h_in_bot, v_h_in);
        const auto in_dihedral_score
            = optimize_dihedral(tri_in_top, tri_in_new) + optimize_dihedral(tri_in_new, tri_in_bot);

        const auto sharpness_tri_in_top = optimize_angle(to_insert_position, v_h_in, v_h_in_top, opts);
        const auto sharpness_tri_in_new = optimize_angle(to_insert_position, center_position, v_h_in, opts);
        const auto sharpness_tri_in_bot = optimize_angle(center_position, v_h_in_bot, v_h_in, opts);
        const auto in_sharpness = sharpness_tri_in_top + sharpness_tri_in_new + sharpness_tri_in_bot;

        const auto total_out = out_dihedral_score + out_sharpness;
        const auto total_in = in_dihedral_score + in_sharpness;

        h_out_candidates.emplace_back(Candidate{h.halfedge(), tri_out_new, total_out, v_h_out});
        h_in_candidates.emplace_back(Candidate{h.halfedge(), tri_in_new, total_in, v_h_in});
    });

    std::sort(h_out_candidates.begin(), h_out_candidates.end());
    std::sort(h_in_candidates.begin(), h_in_candidates.end());

    // probably don't need to check all 20 - 42 pairs
    for (const auto& h_in : h_in_candidates) {
        for (const auto& h_out : h_out_candidates) {
            if (h_in.h == h_out.h) { continue; }
            const auto score_mid = optimize_dihedral(h_in.normal, h_out.normal);
            const auto valency_score = optimize_valency(m, h_out.h, h_in.h, opts);
            h_pair_candidates.emplace(CandidatePair{
                h_in.h, h_out.h, h_in.score + h_out.score + score_mid + valency_score
            });
        }
    }
    edges[0] = h_pair_candidates.top().out;
    edges[1] = m.walker(h_pair_candidates.top().in).opp().halfedge();

    return edges;
}


inline auto reexpand_points(::HMesh::Manifold& manifold, Collapse&& collapse, const std::vector<Vec3>& points,
                            const ReexpandOptions& opts = ReexpandOptions()) -> void
{
    std::cout << "reexpanding\n";
    const auto& manifold_positions = manifold.positions;

    std::unordered_multimap<Point, VertexID, PointHash, PointEquals> point_to_manifold_ids;
    for (auto manifold_vid : manifold.vertices()) {
        auto pos = manifold_positions[manifold_vid];
        point_to_manifold_ids.emplace(pos, manifold_vid);
    }

    auto real_index_to_manifold_iter = [&](const NodeID id) {
        const auto pos = points[id];
        auto [fst, snd] = point_to_manifold_ids.equal_range(pos);
        return std::ranges::subrange(fst, snd) | std::views::values;
    };
    using IndexIter = decltype(real_index_to_manifold_iter(std::declval<NodeID>()));

    auto try_two_edge_expand = [&](IndexIter vs, const Point& p) -> bool {
        // we want to get as close to 90 degrees as possible here
        for (auto this_vert : vs) {
            const auto candidate = find_edge_pair(manifold, this_vert, p, opts);
            if (candidate[0] != InvalidHalfEdgeID && candidate[1] != InvalidHalfEdgeID) {
                const auto vnew = manifold.split_vertex(candidate[1], candidate[0]);
                GEL_ASSERT_NEQ(vnew, InvalidVertexID);
                manifold.positions[vnew] = p;
                point_to_manifold_ids.emplace(p, vnew);
                return true;
            }
        }
        return false;
    };

    size_t failures = 0;
    for (auto collapse_iter : collapse | std::views::reverse) {
        for (auto [active, latent] : collapse_iter | std::views::reverse) {
            // find the manifold_ids for the active vertex
            const auto manifold_ids = real_index_to_manifold_iter(active);
            const auto point_position = points[latent];

            if (try_two_edge_expand(manifold_ids, point_position)) continue;

            failures++;
        }
    }
    std::cerr << "failures: " << failures << "\n";
}
} // namespace HMesh

#endif // GEL_HMESH_COLLAPSE_H
