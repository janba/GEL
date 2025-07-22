//
// Created by Cem Akarsubasi on 7/11/25.
// Some Collapse data structures and functions that can be separated from RsR

#ifndef GEL_HMESH_COLLAPSE_H
#define GEL_HMESH_COLLAPSE_H

#include <GEL/HMesh/Manifold.h>
#include <GEL/Util/Assert.h>
#include <GEL/Util/ParallelAdapters.h>
#include <span>
#include <unordered_map>
#include <vector>
#include <numbers>

namespace HMesh {

    using Vec3 = CGLA::Vec3d;
    using Point = Vec3;
    using NodeID = size_t;

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
        friend struct CollapseTest;
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
            std::unordered_map<ActiveNodeID, LatentNodeID> activity;

        public:
            auto insert(ActiveNodeID active, LatentNodeID latent) -> void
            {
                GEL_ASSERT_NEQ(active, latent);
                GEL_ASSERT_NEQ(active, InvalidNodeID);
                GEL_ASSERT_NEQ(latent, InvalidNodeID);
                GEL_ASSERT(!activity.contains(active));
                GEL_ASSERT(!activity.contains(latent));
                activity.emplace(active, latent);
                activity.emplace(latent, Collapse::InvalidNodeID);
            }

            auto is_used(const NodeID id) const -> bool { return activity.contains(id); }

            auto is_active(const ActiveNodeID id) const -> bool
            {
                if (const auto maybe = activity.find(id); maybe != activity.end()) {
                    return maybe->second != Collapse::InvalidNodeID;
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
            for (auto [active, latent]: activity_map.activity) {
                if (latent != InvalidNodeID) {
                    m_collapses.emplace_back(active, latent);
                    ++items;
                }
            }
            m_collapse_ranges.emplace_back(CollapseInfo{begin_idx, begin_idx + items});

            // update remaining indices
            auto [fst, snd] = std::ranges::remove_if(m_remaining, [&](const NodeID id) {
                const auto maybe = activity_map.activity.find(id);
                return (maybe != activity_map.activity.end() && maybe->second == Collapse::InvalidNodeID);
            });
            m_remaining.erase(fst, snd);
        }

        [[nodiscard]]
        auto get_collapse_span(const size_t at) const -> std::span<const std::pair<NodeID, NodeID>>
        {
            if (at >= m_collapse_ranges.size()) { throw std::out_of_range("collapse_ranges"); }
            return {m_collapses.begin() + m_collapse_ranges.at(at).begin,
                    m_collapses.begin() + m_collapse_ranges.at(at).end};
        }
    };
    static_assert(std::ranges::viewable_range<Collapse&>);
    static_assert(std::ranges::viewable_range<const Collapse&>);

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

    inline CGLA::Vec3d half_edge_direction(const HMesh::Manifold& m, HMesh::HalfEdgeID h)
    {
        const auto w = m.walker(h);
        const auto current = w.vertex();
        const auto opposing = w.opp().vertex();
        return CGLA::normalize(m.positions[opposing] - m.positions[current]);
    }

    inline Vec3d triangle_normal(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3)
    {
        const auto v1 = p2 - p1;
        const auto v2 = p3 - p1;
        return CGLA::normalize(CGLA::cross(v1, v2));
    }

    // returns 0 at 180 degrees, 1 at 90 (or 270) degrees
    inline double optimize_dihedral(const Vec3d n1, const Vec3d n2)
    {
        const auto angle = CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0;
        return std::abs(angle);
    }

    // returns 0 for an equilateral triangle,
    inline double optimize_angle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3, const ReexpandOptions& opts)
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
        constexpr auto inf = std::numeric_limits<double>::infinity();
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
                                                           const CGLA::Vec3d& to_insert_position, const ReexpandOptions& opts)
    {
        std::array edges = {HMesh::InvalidHalfEdgeID, HMesh::InvalidHalfEdgeID};
        const auto center_position = m.positions[center_idx];

        struct Candidate {
            HalfEdgeID h;
            Vec3d normal;
            double score;
            Vec3d points;
            std::weak_ordering operator<=>(const Candidate& other) const
            {
                return this->score == other.score  ? std::weak_ordering::equivalent
                       : this->score < other.score ? std::weak_ordering::less
                                                   : std::weak_ordering::greater;
            }
        };
        struct CandidatePair {
            HalfEdgeID in;
            HalfEdgeID out;
            double score;
            std::weak_ordering operator<=>(const CandidatePair& other) const
            {
                return this->score == other.score  ? std::weak_ordering::equivalent
                       : this->score < other.score ? std::weak_ordering::less
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
        for (const auto& h_in: h_in_candidates) {
            for (const auto& h_out: h_out_candidates) {
                if (h_in.h == h_out.h) { continue; }
                const auto score_mid = optimize_dihedral(h_in.normal, h_out.normal);
                const auto valency_score = optimize_valency(m, h_out.h, h_in.h, opts);
                h_pair_candidates.emplace(CandidatePair{h_in.h, h_out.h, h_in.score + h_out.score + score_mid + valency_score});
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
        for (auto manifold_vid: manifold.vertices()) {
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
            for (auto this_vert: vs) {
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
        for (auto collapse_iter: collapse | std::views::reverse) {
            for (auto [active, latent]: collapse_iter | std::views::reverse) {
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
