//
// Created by Cem Akarsubasi on 9/10/25.
//

#include "HierarchicalReconstruction.h"

#include <unordered_map>

namespace HMesh::RSR
{

using namespace detail;
using namespace Util::detail;

bool vec3_eq(const Vec3& lhs, const Vec3& rhs, double eps = 1e-4)
{
    return lhs.all_le(rhs + Vec3(eps)) && lhs.all_ge(rhs - Vec3(eps));
}

bool float_eq(double lhs, double rhs, double eps = 1e-4)
{
    return std::abs(lhs - rhs) < eps;
}

Vec3 half_edge_direction(const HMesh::Manifold& m, HMesh::HalfEdgeID h)
{
    const auto w = m.walker(h);
    const auto current = w.vertex();
    const auto opposing = w.opp().vertex();
    return CGLA::normalize(m.positions[opposing] - m.positions[current]);
}

Vec3 triangle_normal(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    const auto v1 = p2 - p1;
    const auto v2 = p3 - p1;
    return CGLA::normalize(CGLA::cross(v1, v2));
}

// returns 0 at 180 degrees, 1 at 90 (or 270) degrees
double optimize_dihedral(const Vec3& n1, const Vec3& n2)
{
    const auto angle = CGLA::dot(n1, n2) - 1.0;
    return std::abs(angle);
}

// returns 0 for an equilateral triangle
double optimize_angle(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const double& angle_factor,
    const double& angle_threshold_penalty,
    const double& angle_threshold_cos
    )
{
    const auto e1 = p2 - p1;
    const auto e1_len = e1.length();
    const auto e2 = p3 - p1;
    const auto e2_len = e2.length();

    const auto e4 = p3 - p2;
    const auto e4_len   = e4.length();

    std::array angles {
        dot(e1, e2) / (e1_len * e2_len),
        -dot(e2, e4) / (e2_len * e4_len),
        dot(-e4, -e1) / (e4_len * e1_len)
    };

    std::array lengths {e4_len, e1_len, e2_len};
    const auto shortest = std::min(e4_len, std::min(e1_len, e2_len));
    double min_angle_ = 0;
    double opposing_distance = 0;
    for (auto i = 0; i < 3; ++i) {
        if (angles[i] < 0) {
            continue;
        } else if (angles[i] > min_angle_) {
            min_angle_ = angles[i];
            opposing_distance = lengths[i];
        }
    }
    // TODO: this is slow
    auto min_angle_acos = std::acos(min_angle_);
        //1 - min_angle_; // ol reliable
    if (min_angle_ > angle_threshold_cos) {
        return angle_threshold_penalty;
    } else {
        auto score = std::abs(
            std::numbers::pi / 3.0
            //0.5
            - min_angle_acos);
        return score * shortest * angle_factor;
    }
}


struct Split {
    HMesh::HalfEdgeID h_in = InvalidHalfEdgeID;
    HMesh::HalfEdgeID h_out = InvalidHalfEdgeID;
    double max_angle = 0;
};

struct Triangle {
    double shared_length = 0;
    Vec3 normal = {0, 0, 0};

    Triangle() = default;

    Triangle(const Vec3& p1, const Vec3& p2, const Vec3& p3, double shared_length) :
        shared_length(shared_length), normal(triangle_normal(p1, p2, p3))
    {}
};

auto one_ring_max_angle(const Manifold& manifold, const VertexID vid) -> double
{

    auto dihedral_from_hedge = [&manifold](HalfEdgeID h) {
        Walker w = manifold.walker(h);
        auto v1 = w.vertex();
        auto v2 = w.next().vertex();
        auto v3 = w.next().next().vertex();
        auto v4 = w.opp().vertex();
        auto v5 = w.opp().next().vertex();
        auto v6 = w.opp().next().next().vertex();
        if (v1 == InvalidVertexID ||
            v2 == InvalidVertexID ||
            v3 == InvalidVertexID ||
            v4 == InvalidVertexID ||
            v5 == InvalidVertexID ||
            v6 == InvalidVertexID) {
            GEL_ASSERT(false, "wat");
            return 0.0;
        }
        auto n0 = triangle_normal(
            manifold.positions[v1],
            manifold.positions[v2],
            manifold.positions[v3]);
        auto n1 = triangle_normal(
            manifold.positions[v4],
            manifold.positions[v5],
            manifold.positions[v6]);
        return optimize_dihedral(n0, n1);
    };

    double max_angle = 0;
    for (const auto h: manifold.incident_halfedges(vid)) {
        GEL_ASSERT_NEQ(h, InvalidHalfEdgeID);
        max_angle = std::max(max_angle, dihedral_from_hedge(h));
    }
    return max_angle;
}

Split find_edge_pair(const Manifold& m, const VertexID center_idx, const Vec3& v_new_position,
                     const Vec3& v_old_position, const ReexpandOpts& opts, double angle_threshold_cos)
{
    const auto angle_factor = opts.angle_factor;
    const auto angle_threshold_penalty = opts.angle_threshold_penalty;

    // this is getting too complicated
    const auto print_hedge = [&m](HalfEdgeID he) {
        auto v_from = m.positions[m.walker(he).opp().vertex()];
        auto v_to = m.positions[m.walker(he).vertex()];
        std::cout << v_from << " -> " << v_to << " (" << (v_to - v_from) << ")\n";
    };

    // Optimize the dihedral value
    const auto dihedral_one_ring = [](const Triangle& t1, const Triangle& t2) -> std::pair<double, double> {
        auto d = optimize_dihedral(t1.normal, t2.normal);
        auto edge_length = t1.shared_length;
        return std::make_pair(d, edge_length);
    };

    struct ExpandResult {
        double score = INFINITY;
        double max_angle = INFINITY;
    };
    const auto expand_score = [&](
        const HalfEdgeID& h_in_opp,
        const HalfEdgeID& h_out,
        const Point& v_new_position,
        const Point& v_old_position) -> ExpandResult {
        const auto walker_out = m.walker(h_out);
        const auto walker_in_opp = m.walker(h_in_opp);
        const auto& v_h_out = m.pos(walker_out.vertex());
        const auto& v_h_in  = m.pos(walker_in_opp.vertex());

        constexpr double EPS = 1e-8;

        const double  in_len = std::max((v_old_position - v_h_in).length(), EPS);
        const double out_len = std::max((v_new_position - v_h_out).length(), EPS);
        const auto tri_center_in = Triangle(v_old_position, v_new_position, v_h_in, in_len);
        const auto tri_center_out = Triangle(v_new_position, v_old_position, v_h_out, out_len);

        double angle_cost = 0;
        if (angle_factor > 0.0) {
            angle_cost += optimize_angle(v_old_position, v_new_position, v_h_in, angle_factor, angle_threshold_penalty, angle_threshold_cos);
            angle_cost += optimize_angle(v_new_position, v_old_position, v_h_out,angle_factor, angle_threshold_penalty, angle_threshold_cos);
        }

        double dihedral_cost = 0;
        double max_angle = 0;
        std::array<Triangle, 2> triangles_buffer_;
        auto calculate_one_dihedral = [&](const Triangle& tri) {
            triangles_buffer_[0] = triangles_buffer_[1];
            triangles_buffer_[1] = tri;
            auto [angle, edge_length] = dihedral_one_ring(triangles_buffer_[0], triangles_buffer_[1]);
            dihedral_cost += angle * angle * edge_length;
            max_angle = std::max(max_angle, angle);
        };

        auto walker_prev = walker_out;
        auto walker_next = walker_prev.circulate_vertex_ccw();
        // from out towards in counterclockwise
        triangles_buffer_[1] = tri_center_out;
        while (walker_prev.halfedge() != walker_in_opp.halfedge()) {
            const auto& p2 = m.pos(walker_prev.vertex());
            const auto& p3 = m.pos(walker_next.vertex());
            const auto shared_length = std::max((v_new_position - p3).length(), EPS);
            Triangle tri = {v_new_position, p2, p3, shared_length};
            if (angle_factor > 0.0) {
                angle_cost += optimize_angle(v_new_position, p2, p3, angle_factor, angle_threshold_penalty, angle_threshold_cos);
            }
            calculate_one_dihedral(tri);

            walker_prev = walker_prev.circulate_vertex_ccw();
            walker_next = walker_next.circulate_vertex_ccw();
        }
        calculate_one_dihedral(tri_center_in);

        // from in towards out counterclockwise
        walker_prev = walker_in_opp;
        walker_next = walker_in_opp.circulate_vertex_ccw();
        while (walker_prev.halfedge() != walker_out.halfedge()) {
            const auto& p2 = m.pos(walker_prev.vertex());
            const auto& p3 = m.pos(walker_next.vertex());
            const auto shared_length = std::max((v_old_position - p3).length(), EPS);
            Triangle tri = {v_old_position, p2, p3, shared_length};
            if (angle_factor > 0.0) {
                angle_cost += optimize_angle(v_old_position, p2, p3, angle_factor, angle_threshold_penalty, angle_threshold_cos);
            }
            calculate_one_dihedral(tri);

            walker_prev = walker_prev.circulate_vertex_ccw();
            walker_next = walker_next.circulate_vertex_ccw();
        }
        calculate_one_dihedral(tri_center_out);

        const auto [dihedral0_angle, dihedral0_length] = dihedral_one_ring(tri_center_in, tri_center_out);
        const auto dihedral0 = dihedral0_angle * dihedral0_angle * dihedral0_length;
        max_angle = std::max(max_angle, dihedral0_angle);

        return ExpandResult{
            .score = dihedral_cost + dihedral0 + angle_cost,
            .max_angle = max_angle,
        };
    };

    double min_score = INFINITY;
    double max_angle = INFINITY;
    HalfEdgeID h_in_opp;
    HalfEdgeID h_out;
    constexpr double dihedral_threshold = 1.75;

    double min_score_alt = INFINITY;
    double max_angle_alt = INFINITY;
    HalfEdgeID h_in_opp_alt = InvalidHalfEdgeID;
    HalfEdgeID h_out_alt = InvalidHalfEdgeID;

    // TODO: ruling out "bad pairs" early can speed this up by a lot
    auto v_bar = m.positions[center_idx];
    auto ref_length = (v_old_position - v_new_position).length();
    auto norm = m.normal(center_idx);
    auto v3 = v_bar + norm * ref_length;

    for (auto h1 : m.incident_halfedges(center_idx)) {
        auto v1 = m.positions[m.walker(h1).vertex()];
        for (auto h2 : m.incident_halfedges(center_idx)) {
            // Illegal
            if (h1 == h2) {
                continue;
            }
            auto v2 = m.positions[m.walker(h2).vertex()];
            auto plane_norm = CGLA::normalize(CGLA::cross(v3 - v1, v2 - v1));
            //if (CGLA::dot((v_new_position - v_bar), plane_norm) < 0.0 || CGLA::dot((v_old_position - v_bar), plane_norm) > 0.0)
            //    continue;


            auto score = expand_score(h1, h2, v_new_position, v_old_position);
            if (score.score < min_score) {
                min_score = score.score;
                max_angle = score.max_angle;
                h_in_opp = h1;
                h_out = h2;
            }
            if (score.score < min_score_alt && score.max_angle < dihedral_threshold) {
                min_score_alt = score.score;
                max_angle_alt = score.max_angle;
                h_in_opp_alt = h1;
                h_out_alt = h2;
            }
            if (opts.debug_mask & RE_SPLITS) {
                std::cout << "h1: ";
                print_hedge(h1);
                std::cout << "h2: ";
                print_hedge(h2);
                std::cout << "score    : " << score.score << "\n";
                std::cout << "max angle: " << score.max_angle << "\n";
            }
        }
    }
    if (max_angle > dihedral_threshold && max_angle_alt < dihedral_threshold) {
        if (opts.debug_mask & RE_SPLIT_RESULTS) {
            std::cout << "using alternative split\n";
            std::cout << "h_in_opp: "; print_hedge(h_in_opp_alt);
            std::cout << "h_out:    "; print_hedge(h_out_alt);

        }
        if (h_out_alt == InvalidHalfEdgeID) {
            return {};
        }
        return Split{m.walker(h_in_opp_alt).opp().halfedge(), h_out_alt, max_angle_alt};
    }
    if (opts.debug_mask & RE_SPLIT_RESULTS) {
        std::cout << "using regular split\n";
        std::cout << "h_in_opp: "; print_hedge(h_in_opp);
        std::cout << "h_out:    "; print_hedge(h_out);
    }
    if (h_out == InvalidHalfEdgeID) {
        return {};
    }
    return Split{m.walker(h_in_opp).opp().halfedge(), h_out, max_angle};
}

auto angle_flip_check(const Manifold& manifold, const HalfEdgeID he, double angle_threshold) -> bool
{
    if (!manifold.precond_flip_edge(he))
        return false;
    const auto angles = [&](HalfEdgeID he) -> std::pair<double, double> {
        auto walker = manifold.walker(he);
        auto v1 = manifold.positions[walker.vertex()];
        auto v2 = manifold.positions[walker.next().vertex()];
        auto v3 = manifold.positions[walker.opp().vertex()];
        auto e_shared = CGLA::normalize(v1 - v3);
        auto e_next = CGLA::normalize(v2 - v1);
        auto e_prev = CGLA::normalize(v2 - v3);

        auto angle1 = CGLA::dot(e_next, -e_shared);
        auto angle2 = CGLA::dot(e_prev, e_shared);

        return {angle1, angle2};
    };
    const auto he_opp = manifold.walker(he).opp().halfedge();
    auto [angle1, angle2] = angles(he);
    auto [angle3, angle4] = angles(he_opp);

    if (angle1 < 0.0 || angle2 < 0.0 || angle3 < 0.0 || angle4 < 0.0) {
        return false;
    }

    if (std::acos(angle1) < angle_threshold ||
        std::acos(angle2) < angle_threshold ||
        std::acos(angle3) < angle_threshold ||
        std::acos(angle4) < angle_threshold) {
        return true;
    }
    return false;
}

auto collapse_points(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
                     const CollapseOpts& options) -> std::pair<Collapse, PointCloud>
{
    if (options.max_iterations == 0) {
        return std::make_pair(Collapse(), PointCloud(vertices, normals));
    }
    std::cout << "Collapsing..." << std::endl;
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    ImmediatePool pool;
    QuadraticCollapseGraph graph;

    // initialize graph
    for (auto i = 0UL; i < vertices.size(); ++i) {
        graph.add_node(vertices[i], normals[i]);
    }
    auto indices = [&vertices] {
        std::vector<NodeID> temp(vertices.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();

    const auto kd_tree = build_kd_tree_of_indices(vertices, indices);
    const auto neighbor_map = calculate_neighbors(pool, vertices, kd_tree, options.initial_neighbors);

    // This also initializes distances
    for (const auto& neighbors : neighbor_map) {
        const NodeID this_id = neighbors[0].id;
        for (const auto& neighbor : neighbors | std::views::drop(1)) {
            // kNN connection
            graph.connect_nodes(this_id, neighbor.id);
        }
    }

    std::vector<std::vector<SingleCollapse>> collapses;

    for (size_t iter = 0; iter < options.max_iterations; ++iter) {
        // TODO: stricter checking
        const size_t max_collapses =
            [&]() -> size_t {
                if (options.max_collapses != 0) return options.max_collapses;

                return vertices.size() * std::pow(0.5, iter) * options.reduction_per_iteration;
            }();

        std::vector<SingleCollapse> activity;

        size_t count = 0;
        while (count < max_collapses) {
            count++;
            auto [active, latent, active_point_coords, latent_point_coords, v_bar] = graph.collapse_one();

            activity.emplace_back(active_point_coords, latent_point_coords, v_bar);
        }
        collapses.emplace_back(std::move(activity));
        std::cout << "Collapsed " << count << " of " << max_collapses << std::endl;
    }
    return std::make_pair(Collapse {std::move(collapses)}, graph.to_point_cloud()); // TODO
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

std::optional<HalfEdgeID> find_crossed_edge(
    const Manifold& manifold,
    const VertexID id,
    const Point& active_pos,
    const Point& latent_pos,
    const ReexpandOpts& opts)
{
    for (auto he: manifold.incident_halfedges(id)) {
        auto walker = manifold.walker(he);
        auto flip_maybe = walker.next().halfedge();
        auto v1 = manifold.positions[walker.vertex()];
        auto v2 = manifold.positions[walker.next().vertex()];
        auto tangent = normalize(v1 - v2);
        if (opts.debug_mask & RE_FIRST_FLIP) {
            std::cout << v1 << "\n";
            std::cout << v2 << "\n";
        }

        auto normal = CGLA::normalize(CGLA::cross((active_pos - v1), tangent));

        // normal of the plane
        Vec3 binormal = CGLA::cross(tangent, normal);
        // point on the plane
        const Point& p0 = (v1 + v2) * 0.5;
        const double r = (v1 - v2).length() * 0.5;

        // starting point in our line
        const Point& l0 = active_pos;
        const Vec3 l = normalize(latent_pos - active_pos);

        auto den = dot(l, binormal);
        if (std::abs(den) < 1e-8) {
            continue;
        } // parallel case
        auto len = (latent_pos - active_pos).length();
        auto d = dot((p0 - l0), binormal) / den;

        auto p = l0 + d * l;
        auto distance = (p - p0).length();
        if (opts.debug_mask & RE_FIRST_FLIP) {
            std::cout << "d  : " << d << "\n";
            std::cout << "len: " << len << "\n";
        }
        if (distance > r) {
            continue;
        }

        if (d > 0 && len > (d * 0.90) && manifold.precond_flip_edge(flip_maybe)) {
            if (opts.debug_mask & RE_FIRST_FLIP) {
                std::cout << "flipped!" << "\n";
            }
            return flip_maybe;
        }
    }
    return std::nullopt;
}

// Fast quick sort implementations from cppreference
namespace
{
    template <typename Comp = std::less<>>
    void quick_sort(std::forward_iterator auto first, std::forward_iterator auto last, Comp comp = std::less<>{})
    {
        if (first == last)
            return;

        auto pivot = *std::next(first, std::distance(first, last) / 2);
        auto middle1 = std::partition(first, last, [&](const auto& elem)
        {
            return comp(elem, pivot);
        });
        auto middle2 = std::partition(middle1, last, [&](const auto& elem)
        {
            return !comp(pivot, elem);
        });

        quick_sort(first, middle1, comp);
        quick_sort(middle2, last, comp);
    }

    template <typename Range, typename Comp = decltype(std::less<>{})>
    void quick_sort(Range&& arr, Comp comp = std::less<>{})
    {
        if (std::ranges::distance(arr) == 0) {
            return;
        }
        quick_sort(arr.begin(), arr.end(), comp);
    }
}

auto reexpand_points(Manifold& manifold, Collapse&& collapse, const ReexpandOpts& opts) -> void
{
    std::cout << "reexpanding" << std::endl;
    const auto& manifold_positions = manifold.positions;

    // TODO: Replace this with the collection in Util
    // TODO: factor this out
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

    double angle_threshold_cos = std::cos(opts.angle_threshold);
    // Now we need to consider two position candidates

    size_t expansion_failures = 0;
    size_t bad_expansions = 0;
    size_t flips = 0;
    int iteration = 0;
    std::vector<HalfEdgeID> one_ring;
    std::vector<HalfEdgeID> circle;
    std::vector<HalfEdgeID> two_ring;
    for (const auto& collapse_iter : collapse.collapses | std::views::reverse) {
        for (auto single_collapse : collapse_iter | std::views::reverse) {
            iteration++;
            // find the manifold_ids for the active vertex

            const auto active_pos = single_collapse.active_point_coords;
            const auto latent_pos = single_collapse.latent_point_coords;
            const auto v_bar = single_collapse.v_bar;
            // FIXME: Debug info
            if (opts.debug_mask & RE_ITERATION) {
                std::cout << "--------------------------------\n";
                std::cout << "@iteration: " << iteration << "\n";
                std::cout << "active pos: " << active_pos << "\n";
                std::cout << "latent pos: " << latent_pos << "\n";
                std::cout << "combin pos: " << v_bar << "\n";
            }
            if (opts.debug_mask & RE_MARK_SPLITS) {
                manifold.add_face({active_pos, latent_pos, v_bar});
            }
            if (opts.debug_mask & RE_FIRST_FLIP) {
                std::cout << "First flip:\n";
            }
            // repair local geometry maybe
            const auto manifold_ids = position_to_manifold_iter(v_bar);
            for (const auto id : manifold_ids) {
                auto maybe1 = find_crossed_edge(manifold, id, latent_pos, active_pos, opts);
                if (maybe1) {
                    manifold.flip_edge(*maybe1);
                }
                auto maybe2 = find_crossed_edge(manifold, id, active_pos, latent_pos, opts);
                if (maybe2) {
                    manifold.flip_edge(*maybe2);
                }

                manifold.positions[id] = active_pos;
            }

            // actually do the expansion
            const auto new_vid = [&]() -> VertexID {
                // we want to get as close to 90 degrees as possible here
                for (const auto this_vert : manifold_ids) {
                    const auto candidate = find_edge_pair(manifold, this_vert, latent_pos, active_pos, opts, angle_threshold_cos);
                    if (candidate.h_in != InvalidHalfEdgeID) {
                        const auto vnew = manifold.split_vertex(candidate.h_in, candidate.h_out);
                        GEL_ASSERT_NEQ(vnew, InvalidVertexID);
                        return vnew;
                    }
                }
                return InvalidVertexID;
            }();

            if (new_vid == InvalidVertexID) {
                expansion_failures++;
                continue;
            }

            manifold.positions[new_vid] = latent_pos;

            if (std::ranges::distance(manifold_ids) < 4) {
                InplaceVector<VertexID> copy(manifold_ids.begin(), manifold_ids.end());
                for (auto id: copy) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
            } else {
                std::vector copy(manifold_ids.begin(), manifold_ids.end());
                for (auto id: copy) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
            }
            point_to_manifold_ids.emplace(latent_pos, new_vid);
            point_to_manifold_ids.erase(v_bar);

            if (opts.debug_mask & RE_SECOND_FLIP) {
                std::cout << "Second flip:\n";
            }


            one_ring.clear();
            circle.clear();
            two_ring.clear();
            circulate_vertex_ccw(manifold, new_vid, [&](Walker& w) {
                const HalfEdgeID one_ring_he = w.halfedge();
                const HalfEdgeID circle_he = w.next().halfedge();
                const HalfEdgeID two_ring_he = w.next().opp().next().halfedge();
                const HalfEdgeID two_ring_he2 = w.next().opp().prev().halfedge();

                one_ring.push_back(one_ring_he);
                circle.push_back(circle_he);
                two_ring.push_back(two_ring_he);
                two_ring.push_back(two_ring_he2);
            });

            // sort from longest to shortest edge
            const auto& manifold_cref = manifold;
            auto cmp = [&manifold_cref](const HalfEdgeID& e1, const HalfEdgeID& e2) -> bool {
                // TODO: There is a serious soundness issue here
                GEL_ASSERT(manifold_cref.in_use(e1), "%ld", e1.get_index());
                GEL_ASSERT(manifold_cref.in_use(e2), "%ld", e2.get_index());
                auto len1 = manifold_cref.length(e1);
                auto len2 = manifold_cref.length(e2);
                return len1 < len2;
            };

            for (int i = 0; i < opts.refinement_iterations; ++i) {
                auto threshold = opts.refinement_angle_threshold;
                quick_sort(one_ring, cmp);
                for (HalfEdgeID h : one_ring | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, M_PI / 180.0 * threshold);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
                quick_sort(circle, cmp);
                for (HalfEdgeID h : circle | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, M_PI / 180.0 * threshold);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
                quick_sort(two_ring, cmp);
                for (HalfEdgeID h : two_ring | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, M_PI / 180.0 * threshold);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
            }

            if (opts.debug_mask & RE_ERRORS) {
                const auto v_new_max_angle = one_ring_max_angle(manifold, new_vid);
                const auto v_old_max_angle = one_ring_max_angle(manifold, manifold_ids.front());
                if (v_new_max_angle > opts.angle_stop_threshold || v_old_max_angle > opts.angle_stop_threshold) {
                    bad_expansions++;
                    if (opts.debug_mask & RE_ERRORS) {
                        auto normal = manifold.normal(new_vid);
                        manifold.add_face({active_pos, latent_pos, v_bar + (active_pos-latent_pos).length() * normal});
                        std::cout << "v_new max angle: " << v_new_max_angle << "\n";
                        std::cout << "v_old max angle: " << v_old_max_angle << "\n";
                        std::cout << "failed at iteration " << iteration << "\n";
                    }
                    if (opts.early_stop_at_error && opts.stop_at_error == 0) {
                        goto EXIT;
                    }
                    if (opts.stop_at_error > 0 && opts.stop_at_error == bad_expansions) {
                        goto EXIT;
                    }
                }
            }

            if (opts.stop_at_iteration > 0 && iteration == opts.stop_at_iteration) {
                std::cout << "stopped early at " << iteration << "\n";
                goto EXIT;
            }
        }
    }
    EXIT:
    std::cout << "flips: " << flips << "\n";
    std::cerr << "failures: " << expansion_failures << "\n";
}

}
