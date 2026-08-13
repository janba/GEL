//
// Created by Cem Akarsubasi on 9/10/25.
//

#include <GEL/Util/InplaceVector.h>

#include <GEL/HMesh/HierarchicalReconstruction.h>

#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/obj_save.h>
#include <numbers>
#include <unordered_map>

namespace HMesh::RSR
{
using NodeID = size_t;
using Point = CGLA::Vec3d;
using Vec3 = CGLA::Vec3d;
using Vec2 = CGLA::Vec2d;

using namespace detail;
using namespace Util::detail;

using Geometry::AMGraph;

/// Raw collapse is used by the CollapseGraph which
/// needs to track the NodeIDs in addition to the coordinates.
struct RawCollapse {
    /// The point whose ID remains in the graph after the collapse
    NodeID active;
    /// The point whose ID is removed from the graph after the collapse
    NodeID latent;
    /// The old coordinates of active
    Point active_point_coords;
    /// The old coordinates of latent
    Point latent_point_coords;
    /// The new coordinates of the collapsed point.
    Point v_bar;
};

/// Derives from AMGraph to perform the collapse. The actual edges are
/// stored inside a BTree which serves as a mutable priority queue.
struct CollapseGraph : AMGraph {
private:
    struct Vertex {
        Point position;
        Vec3 normal;
        /// The weight of a vertex depends on how many vertices have been collapsed
        /// into it. This penalizes continuously collapsing into the same vertex
        /// and presumably creates a more balanced graph.
        double weight = 1;
    };

    /// Stores edge information
    struct Edge {
        NodeID from;
        NodeID to;
        double dist;

        bool operator==(const Edge& rhs) const
        {
            return (from == rhs.from && to == rhs.to) || (from == rhs.to && to == rhs.from);
        }

        bool operator<(const Edge& other) const
        {
            return dist < other.dist;
        }
    };

public:
    /// Stores vertex information
    Util::AttribVec<NodeID, Vertex> m_vertices;

private:
    /// The collapse queue stored as a BTreeSet
    Util::BTreeSet<Edge> m_collapse_queue;
    /// Map from an EdgeID to the edge length. Since edges are stored sorted
    /// in terms of edge length in m_collapse_queue, extracting them in logarithmic
    /// time requires us to store the edge length in two separate places.
    Util::AttribVec<EdgeID, double> m_edges;

public:
    /// Insert a vertex
    auto add_node(const Point& position, const Vec3& normal) -> NodeID
    {
        const NodeID n = AMGraph::add_node();
        m_vertices[n] = Vertex{.position = position, .normal = normal};
        return n;
    }

    /// Insert an edge into the graph and the priority queue
    auto connect_nodes(const NodeID n0, const NodeID n1) -> EdgeID
    {
        if (n1 > n0) {
            const EdgeID e = AMGraph::connect_nodes(n0, n1);
            auto [_, dist] = distance_function(n0, n1);
            m_collapse_queue.emplace(n0, n1, dist);
            m_edges[e] = dist;
            return e;
        }
        return InvalidEdgeID;
    }

    /// Perform a single collapse. The edge with the lowest "weight" is removed from
    /// the graph and the ends of the edge are combined into a single vertex. All of
    /// the edges that used to belong to the two vertices are removed from the priority
    /// queue and are reinserted with the new coordinates and vertex weights.
    auto collapse_one() -> RawCollapse
    {
        while (!m_collapse_queue.empty()) {
            auto edge_ = m_collapse_queue.begin();
            const auto edge = *edge_;
            m_collapse_queue.erase(edge_);
            const auto active = edge.from;
            const auto latent = edge.to;

            if (AMGraph::find_edge(active, latent) != AMGraph::InvalidEdgeID) {
                const auto active_coords = m_vertices[active].position;
                const auto latent_coords = m_vertices[latent].position;
                GEL_ASSERT_FALSE(active_coords.any([](auto d){ return std::isnan(d); }));
                GEL_ASSERT_FALSE(latent_coords.any([](auto d){ return std::isnan(d); }));

                // recalculate current edges
                for (auto v : AMGraph::neighbors_lazy(active)) {
                    auto v0 = std::min(v, active);
                    auto v1 = std::max(v, active);
                    auto edge_id = find_edge(v0, v1);
                    auto dist = m_edges[edge_id];
                    m_collapse_queue.extract(Edge{v0, v1, dist});
                }

                for (auto v : AMGraph::neighbors_lazy(latent)) {
                    auto v0 = std::min(v, latent);
                    auto v1 = std::max(v, latent);
                    auto edge_id = find_edge(v0, v1);
                    auto dist = m_edges[edge_id];
                    m_collapse_queue.extract(Edge{v0, v1, dist});
                }

                // combine the two vertices
                double total_weight = m_vertices[latent].weight + m_vertices[active].weight;
                const auto new_normal =
                    lerp(m_vertices[latent].normal, m_vertices[active].normal,
                         m_vertices[active].weight / total_weight);
                const auto v_bar =
                    lerp(m_vertices[latent].position, m_vertices[active].position,
                         m_vertices[active].weight / total_weight);
                m_vertices[active].position = v_bar;
                m_vertices[active].normal = new_normal;
                m_vertices[active].weight += m_vertices[latent].weight;
                // we set the latent positions to NaN to defend against reusing them.
                m_vertices[latent].position = Point(std::numeric_limits<double>::signaling_NaN());
                m_vertices[latent].normal = Vec3(std::numeric_limits<double>::signaling_NaN());

                // recalculate current edges
                for (auto v : AMGraph::neighbors_lazy(active)) {
                    auto [v0, v1] = std::minmax(v, active);
                    connect_nodes(v0, v1);
                }

                for (auto v : AMGraph::neighbors_lazy(latent)) {
                    auto [v0, v1] = std::minmax(v, active);
                    connect_nodes(v0, v1);
                }
                AMGraph::erase_node(latent);

                return RawCollapse{
                    .active = active,
                    .latent = latent,
                    .active_point_coords = active_coords,
                    .latent_point_coords = latent_coords,
                    .v_bar = v_bar
                };
            }
        }
        // Immediately fail if we run out of edges to collapse. This might happen if the graph is initialized
        // with an insufficient number of initial nearest neighbors. Not a hugely important case to handle right now.
        GEL_ASSERT(false, "Ran out of edges to collapse");
        return RawCollapse{};
    }

    /// Return the remaining points
    auto to_point_cloud() -> PointCloud
    {
        std::vector<Point> points;
        std::vector<Vec3> normals;
        for (auto i = 0UL; i < m_vertices.size(); ++i) {
            if (!m_vertices[i].position.any([](const double e) { return std::isnan(e); })) {
                points.emplace_back(m_vertices[i].position);
                normals.emplace_back(m_vertices[i].normal);
            }
        }
        return PointCloud{std::move(points), std::move(normals)};
    }

private:
    /// Returns the optimal point and the collapse distance for two vertices
    [[nodiscard]]
    std::pair<CGLA::Vec3d, double> distance_function(const NodeID n0, const NodeID n1) const
    {
        GEL_ASSERT(valid_node_id(n0) && valid_node_id(n1));
        if (valid_node_id(n0) && valid_node_id(n1)) {
            const auto& v0 = m_vertices[n0];
            const auto& v1 = m_vertices[n1];

            const auto total_weight = v0.weight + v1.weight;
            const auto center = lerp(v0.position, v1.position, v0.weight / total_weight);
            const auto tangent_distance =
                Geometry::tangent_space_distance(v0.position - v1.position, v0.normal, v1.normal);

            return std::make_pair(center, tangent_distance * total_weight);
        } else {
            return std::make_pair(CGLA::Vec3d(), CGLA::CGLA_NAN);
        }
    }
};

/// Normalized normal vector of a triangle
Vec3 triangle_normal(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    const auto v1 = p2 - p1;
    const auto v2 = p3 - p1;
    return CGLA::normalize(CGLA::cross(v1, v2));
}

/// returns 0 at 180 degrees, 1 at 90 (or 270) degrees
double optimize_dihedral_angle(const Vec3& n1, const Vec3& n2)
{
    const auto angle = CGLA::dot(n1, n2) - 1.0;
    return std::abs(angle);
}

/// returns 0 for an equilateral triangle
double optimize_min_angle(
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    double angle_factor,
    double angle_threshold_penalty,
    double angle_threshold_cos
)
{
    const auto e1 = p2 - p1;
    const auto e1_len = e1.length();
    const auto e2 = p3 - p1;
    const auto e2_len = e2.length();

    const auto e4 = p3 - p2;
    const auto e4_len = e4.length();

    std::array angles{
        dot(e1, e2) / (e1_len * e2_len),
        -dot(e2, e4) / (e2_len * e4_len),
        dot(-e4, -e1) / (e4_len * e1_len)
    };

    std::array lengths{e4_len, e1_len, e2_len};
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
    // Calculating acos directly is actually rather slow
    auto min_angle_acos = std::acos(min_angle_);
    //1 - min_angle_; // ol reliable
    // Attempt to avoid recalculating this. Probably better to get rid of this function altogether
    // than to eliminate this.
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

/// Information about a split
struct Split {
    HMesh::HalfEdgeID h_in = InvalidHalfEdgeID;
    HMesh::HalfEdgeID h_out = InvalidHalfEdgeID;
    /// From 0.0 (0 degrees) to 2.0 (180 degrees)
    double max_dihedral_angle = 0.0;
    bool isSameSide = false;
};

/// An abstraction of a triangle to make it simple to perform a fold operation over a circular
/// list of triangles.
struct Triangle {
    /// The length of the edge that this triangle shares with the next triangle in that list
    double shared_length = 0;
    Vec3 normal = {0, 0, 0};

    Triangle() = default;

    Triangle(const Vec3& p1, const Vec3& p2, const Vec3& p3, double shared_length) :
        shared_length(shared_length), normal(triangle_normal(p1, p2, p3))
    {}
};

/// Maximum dihedral angle within a vertex's one ring
auto one_ring_max_dihedral_angle(const Manifold& manifold, const VertexID vid) -> double
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
        return optimize_dihedral_angle(n0, n1);
    };

    double max_angle = 0;
    for (const auto h : manifold.incident_halfedges(vid)) {
        GEL_ASSERT_NEQ(h, InvalidHalfEdgeID);
        max_angle = std::max(max_angle, dihedral_from_hedge(h));
    }
    return max_angle;
}

double clamp_local(double value, double upper_bound, double lower_bound) {
    double output = value;
    if (output > upper_bound)
        output = upper_bound;
    if (output < lower_bound)
        output = lower_bound;
    return output;
}

double cal_radians_3d_local(const Vec3& branch_vec, const Vec3& normal, const Vec3& ref_vec) {
    Vec3 proj_vec = branch_vec - dot(normal, branch_vec) /
        normal.length() * normal;
    if (std::abs(proj_vec.length()) < 1e-8)
        return 0.;

    Vec3 proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    float value = clamp_local(
        dot(proj_vec, proj_ref) / proj_vec.length() /
        proj_ref.length(), 1, -1);
    double radian = std::acos(value);
    if (dot(CGLA::cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;
    return radian;
}


/// Find the best edge pair for
Split find_edge_pair(const Manifold& m, const VertexID center_idx, const Vec3& v_new_position,
                     const Vec3& v_old_position, const ReexpandOpts& opts, double angle_threshold_cos)
{
    const auto angle_factor = opts.min_angle_weight;
    const auto angle_threshold_penalty = opts.min_angle_threshold_penalty;

    // Print debug information about a halfedge
    const auto print_hedge = [&m](HalfEdgeID he) {
        auto v_from = m.positions[m.walker(he).opp().vertex()];
        auto v_to = m.positions[m.walker(he).vertex()];
        std::cout << v_from << " -> " << v_to << " (" << (v_to - v_from) << ")\n";
    };

    // Optimize the dihedral value
    const auto dihedral_one_ring = [](const Triangle& t1, const Triangle& t2) -> std::pair<double, double> {
        auto d = optimize_dihedral_angle(t1.normal, t2.normal);
        auto edge_length = t1.shared_length;
        return std::make_pair(d, edge_length);
    };

    struct ExpandResult {
        double score = INFINITY;
        double max_angle = INFINITY;
    };
    // Calculate how bad a given expansion would be
    const auto expand_score = [&](
        const HalfEdgeID& h_in_opp,
        const HalfEdgeID& h_out,
        const Point& v_new_position,
        const Point& v_old_position) -> ExpandResult {
        const auto walker_out = m.walker(h_out);
        const auto walker_in_opp = m.walker(h_in_opp);
        const auto& v_h_out = m.pos(walker_out.vertex());
        const auto& v_h_in = m.pos(walker_in_opp.vertex());

        constexpr double EPS = 1e-8;

        const double in_len = std::max((v_old_position - v_h_in).length(), EPS);
        const double out_len = std::max((v_new_position - v_h_out).length(), EPS);
        const auto tri_center_in = Triangle(v_old_position, v_new_position, v_h_in, in_len);
        const auto tri_center_out = Triangle(v_new_position, v_old_position, v_h_out, out_len);

        double angle_cost = 0;
        if (angle_factor > 0.0) {
            angle_cost += optimize_min_angle(v_old_position, v_new_position, v_h_in, angle_factor,
                                             angle_threshold_penalty, angle_threshold_cos);
            angle_cost += optimize_min_angle(v_new_position, v_old_position, v_h_out, angle_factor,
                                             angle_threshold_penalty, angle_threshold_cos);
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
                angle_cost += optimize_min_angle(v_new_position, p2, p3, angle_factor, angle_threshold_penalty,
                                                 angle_threshold_cos);
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
                angle_cost += optimize_min_angle(v_old_position, p2, p3, angle_factor, angle_threshold_penalty,
                                                 angle_threshold_cos);
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
    // Threshold for a "degenerate" dihedral angle
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
            // TODO: early bad pair culling
            auto plane_norm = CGLA::normalize(CGLA::cross(v3 - v1, v2 - v1));
            //if (CGLA::dot((v_new_position - v_bar), plane_norm) < 0.0 || CGLA::dot((v_old_position - v_bar), plane_norm) > 0.0)
            //    continue;

            // This duplication performs "thresholding" for bad dihedral angles
            // TODO: might be better to move thresholding up instead of having a lot of duplication here
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
            if (opts.debug_opts.debug_mask & RE_SPLITS) {
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
        if (opts.debug_opts.debug_mask & RE_SPLIT_RESULTS) {
            std::cout << "using alternative split\n";
            std::cout << "h_in_opp: ";
            print_hedge(h_in_opp_alt);
            std::cout << "h_out:    ";
            print_hedge(h_out_alt);
        }
        if (h_out_alt == InvalidHalfEdgeID) {
            return {};
        }
        return Split{m.walker(h_in_opp_alt).opp().halfedge(), h_out_alt, max_angle_alt};
    }
    if (opts.debug_opts.debug_mask & RE_SPLIT_RESULTS) {
        std::cout << "using regular split\n";
        std::cout << "h_in_opp: ";
        print_hedge(h_in_opp);
        std::cout << "h_out:    ";
        print_hedge(h_out);
    }
    if (h_out == InvalidHalfEdgeID) {
        return {};
    }
    return Split{m.walker(h_in_opp).opp().halfedge(), h_out, max_angle};
}


Split find_edge_pair(const Manifold& m, const VertexID center_idx, const Vec3& v_new_position,
    const Vec3& v_old_position) {

    auto v_bar = m.positions[center_idx];
    auto norm = m.normal(center_idx);
    norm /= norm.length();
    Vec3 v_new_vector = v_new_position - v_bar;
    Point v_new_plane = v_new_position - dot(v_new_vector, norm) * norm;

    std::vector<double> radians;
    std::vector<HalfEdgeID> hids;

    double max_difference = 0.;
    HalfEdgeID hmax1, hmax2;
    int idx1, idx2;
    bool isSameSide = true;

    // Calculate radians first
    for (auto h1 : m.incident_halfedges(center_idx)) {
        auto v1 = m.positions[m.walker(h1).vertex()];
        radians.push_back(cal_radians_3d_local(v1 - v_bar, norm, v_new_plane - v_bar));
        hids.push_back(h1);
    }

    // Search the best pair
    for (int i = 0; i < radians.size(); i++) {
        double max_difference_this_v = 0.;
        bool isSameSide_this = true;
        HalfEdgeID hmax_this_v;
        int maxidx;
        for (int j = 0; j < radians.size(); j++) {
            if (i == j)
                continue;
            double diff = std::abs(radians[i] - radians[j]);
            if (diff > M_PI)
                diff = 2 * M_PI - diff;
            
            // Shortcut
            if (isSameSide_this) {
                if ((radians[i] <= M_PI && radians[j] > M_PI)
                    || (radians[j] <= M_PI && radians[i] > M_PI)) {
                    isSameSide_this = false;
                    max_difference_this_v = diff;
                    hmax_this_v = hids[j];
                    maxidx = j;
                }
                else {
                    if (diff > max_difference_this_v) {
                        max_difference_this_v = diff;
                        hmax_this_v = hids[j];
                        maxidx = j;
                    }
                }
            }
            else {
                if ((radians[i] <= M_PI && radians[j] > M_PI)
                    || (radians[j] <= M_PI && radians[i] > M_PI)) {
                    if (diff > max_difference_this_v) {
                        max_difference_this_v = diff;
                        hmax_this_v = hids[j];
                        maxidx = j;
                    }
                }
            }
        }

        // Update maximum
        if (isSameSide) {
            if (!isSameSide_this) {
                isSameSide = false;
                max_difference = max_difference_this_v;
                hmax1 = hids[i];
                hmax2 = hmax_this_v;
                idx1 = i;
                idx2 = maxidx;
            }
            else {
                if (max_difference_this_v > max_difference) {
                    max_difference = max_difference_this_v;
                    hmax1 = hids[i];
                    hmax2 = hmax_this_v;
                    idx1 = i;
                    idx2 = maxidx;
                }
            }
        }
        else {
            if (isSameSide_this) {

            }
            else {
                if (max_difference_this_v > max_difference) {
                    max_difference = max_difference_this_v;
                    hmax1 = hids[i];
                    hmax2 = hmax_this_v;
                    idx1 = i;
                    idx2 = maxidx;
                }
            }
        }
    }

    // Decide which one is h_in, which one is h_out, align with v_new, v_old.
    // Smaller radian should be h_in
    HalfEdgeID hin, hout;
    if (radians[idx1] < radians[idx2]) {
        hin = hmax1;
        hout = hmax2;
    }
    else {
        hin = hmax2;
        hout = hmax1;
    }

    //std::cout << norm << std::endl;

    return Split{ m.walker(hin).opp().halfedge(), hout, max_difference, isSameSide };
}

/// For the refinement, returns true if an edge flip will not degenerate the mesh and has at least
/// one corner where the angle is below our given threshold angle (in radians).
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

    // Really expensive, perhaps use a fast approximation or do everything in the cosine domain
    if (std::acos(angle1) < angle_threshold ||
        std::acos(angle2) < angle_threshold ||
        std::acos(angle3) < angle_threshold ||
        std::acos(angle4) < angle_threshold) {
        return true;
    }
    return false;
}

void triangle_cosines(double a2, double b2, double c2,
    double& cosA, double& cosB, double& cosC)
{
    const double eps = 1e-12;

    a2 = std::max(a2, eps);
    b2 = std::max(b2, eps);
    c2 = std::max(c2, eps);

    cosA = (b2 + c2 - a2) / (2.0 * sqrt(b2 * c2));
    cosB = (a2 + c2 - b2) / (2.0 * sqrt(a2 * c2));
    cosC = (a2 + b2 - c2) / (2.0 * sqrt(a2 * b2));

    cosA = std::clamp(cosA, -1.0, 1.0);
    cosB = std::clamp(cosB, -1.0, 1.0);
    cosC = std::clamp(cosC, -1.0, 1.0);
}

auto angle_flip_check_new(const Manifold& manifold, const HalfEdgeID he) {
    if (!manifold.precond_flip_edge(he))
        return false;
    auto walker = manifold.walker(he);
    auto v1 = manifold.positions[walker.vertex()];
    auto v2 = manifold.positions[walker.next().vertex()];
    auto v3 = manifold.positions[walker.opp().vertex()];
    auto v4 = manifold.positions[walker.opp().next().vertex()];
    

    // before flip
    auto e_shared = v1 - v3;
    auto e_next = v2 - v1;
    auto e_prev = v2 - v3;
    auto e_opp_next = v4 - v1;
    auto e_opp_prev = v4 - v3;
    double cosa, cosb, cosc;
    triangle_cosines(CGLA::dot(e_shared, e_shared), CGLA::dot(e_next, e_next),
        CGLA::dot(e_prev, e_prev), cosa, cosb, cosc);
    double maximum_cos = std::max({cosa, cosb, cosc});
    triangle_cosines(CGLA::dot(e_shared, e_shared), CGLA::dot(e_opp_next, e_opp_next),
        CGLA::dot(e_opp_prev, e_opp_prev), cosa, cosb, cosc);
    maximum_cos = std::max({ cosa, cosb, cosc, maximum_cos });

    //after flip
    e_shared = v2 - v4;
    e_next = v1 - v2;
    e_prev = v1 - v4;
    e_opp_next = v3 - v2;
    e_opp_prev = v3 - v4;

    triangle_cosines(CGLA::dot(e_shared, e_shared), CGLA::dot(e_next,e_next),
        CGLA::dot(e_prev,e_prev), cosa, cosb, cosc);
    double maximum_cos2 = std::max({ cosa, cosb, cosc });
    triangle_cosines(CGLA::dot(e_shared, e_shared), CGLA::dot(e_opp_next, e_opp_next),
        CGLA::dot(e_opp_prev, e_opp_prev), cosa, cosb, cosc);
    maximum_cos2 = std::max({ cosa, cosb, cosc, maximum_cos2 });

    return (maximum_cos2 < maximum_cos);
}

void export_graph(const CollapseGraph& g, const std::string& out_path)
{
    std::ofstream file(out_path);
    // Write vertices
    file << "# List of geometric vertices\n";
    for (int i = 0; i < g.m_vertices.size(); i++) {
        Point this_coords = g.m_vertices[i].position;
        file << "v " << std::to_string(this_coords[0])
            << " " << std::to_string(this_coords[1])
            << " " << std::to_string(this_coords[2]) << "\n";
    }

    // Write lines
    file << "\n# Line elements\n";
    for (auto i : g.node_ids()) {
        for (auto neighbor : g.neighbors(i)) {
            if (neighbor < i) {
                file << "l " << (i + 1) << " " << (neighbor + 1) << "\n";
            }
        }
    }
}


auto collapse_points(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
                     const CollapseOpts& opts) -> std::pair<Collapse, PointCloud>
{
    if (opts.max_iterations == 0) {
        return std::make_pair(Collapse(), PointCloud(vertices, normals));
    }
    std::cout << "Collapsing..." << std::endl;
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    ImmediatePool pool;
    CollapseGraph graph;

    // initialize graph
    for (auto i = 0UL; i < vertices.size(); ++i) {
        graph.add_node(vertices[i], normals[i]);
    }
    auto indices = [&vertices] {
        std::vector<NodeID> temp(vertices.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();

    Geometry::Tree kd_tree;
    Geometry::build_kd_tree_of_indices(vertices, indices, kd_tree);
    const auto neighbor_map = Geometry::calculate_neighbors(pool, vertices, kd_tree, opts.initial_neighbors);

    // This also initializes distances
    for (const auto& neighbors : neighbor_map) {
        const NodeID this_id = neighbors[0].id;
        for (const auto& neighbor : neighbors | std::views::drop(1)) {
            if (CGLA::dot(normals[neighbor.id], normals[this_id]) < 0.)
                continue;
            // kNN connection
            graph.connect_nodes(this_id, neighbor.id);
        }
    }

    //export_graph(graph, "collapse_graph.obj");

    std::vector<std::vector<SingleCollapse>> collapses;
    size_t total_collapses = 0;
    for (size_t iter = 0; iter < opts.max_iterations; ++iter) {
        // TODO: stricter checking
        const size_t max_collapses =
            [&]() -> size_t {
                return vertices.size() * std::pow((1. - opts.reduction_per_iteration), iter) * opts.reduction_per_iteration;
            //return vertices.size() * std::pow(opts.reduction_per_iteration, iter + 1);
            }();

        std::vector<SingleCollapse> activity;

        size_t count = 0;
        while (count < max_collapses) {
            total_collapses++;
            count++;
            auto [active, latent, active_point_coords, latent_point_coords, v_bar] = graph.collapse_one();

            activity.emplace_back(active_point_coords, latent_point_coords, v_bar);
            /*if (total_collapses == max_collapses) {
                break;
            }*/
        }
        collapses.emplace_back(std::move(activity));
        std::cout << "Collapsed " << count << " of " << max_collapses << std::endl;
        /*if (total_collapses == max_collapses) {
            break;
        }*/
    }
    std::cout << "Collapsed " << total_collapses << " edges" << std::endl;
    Collapse collapse(std::move(collapses));
    return std::make_pair(std::move(collapse), graph.to_point_cloud());
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

std::optional<HalfEdgeID> find_crossed_edge_new(
    const Manifold& manifold,
    const VertexID id,
    const Point& starting_pos,
    const Point& end_pos,
    const ReexpandOpts& opts)
{
    const double eps = 1e-8;
    const double eps_dist = 1e-6;

    for (auto he : manifold.incident_halfedges(id)) {
        auto walker = manifold.walker(he);
        auto flip_maybe = walker.next().halfedge();

        Vec3 v1 = manifold.positions[walker.vertex()];
        Vec3 v2 = manifold.positions[walker.next().vertex()];

        if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
            std::cout << v1 << "\n";
            std::cout << v2 << "\n";
        }

        // ---------------------------------------
        // 1. Construct plane {v1, v2, starting_pos}
        // ---------------------------------------
        Vec3 e = v2 - v1;
        Vec3 n = CGLA::cross(e, starting_pos - v1);
        double n_len = length(n);

        if (n_len < eps) {
            // Degenerate: starting_pos lies on edge line
            continue;
        }
        n /= n_len;

        // ---------------------------------------
        // 2. Project segment onto plane
        // ---------------------------------------
        auto project = [&](const Vec3& p) {
            return p - dot(p - v1, n) * n;
            };

        Vec3 p0 = project(starting_pos);
        Vec3 p1 = project(end_pos);
        Vec3 q0 = v1;
        Vec3 q1 = v2;

        // ---------------------------------------
        // 3. Build 2D frame
        // ---------------------------------------
        Vec3 u = normalize(q1 - q0);
        Vec3 v = normalize(cross(n, u));

        auto to2D = [&](const Vec3& p) {
            Vec3 d = p - q0;
            return Vec2(dot(d, u), dot(d, v));
            };

        Vec2 P0 = to2D(p0);
        Vec2 P1 = to2D(p1);
        Vec2 Q0(0.0, 0.0);
        Vec2 Q1(length(q1 - q0), 0.0);

        // ---------------------------------------
        // 4. 2D intersection test
        // ---------------------------------------
        auto orient = [](const Vec2& a, const Vec2& b, const Vec2& c) {
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            };

        auto on_segment = [&](const Vec2& a, const Vec2& b, const Vec2& p) {
            return std::min(a[0], b[0]) - eps <= p[0] && p[0] <= std::max(a[0], b[0]) + eps &&
                std::min(a[1], b[1]) - eps <= p[1] && p[1] <= std::max(a[1], b[1]) + eps;
            };

        double o1 = orient(P0, P1, Q0);
        double o2 = orient(P0, P1, Q1);
        double o3 = orient(Q0, Q1, P0);
        double o4 = orient(Q0, Q1, P1);

        bool intersect = false;

        // Proper intersection
        if ((o1 * o2 < -eps) && (o3 * o4 < -eps)) {
            intersect = true;
        }

        // Endpoint / collinear cases
        if (!intersect) {
            if (std::abs(o1) < eps && on_segment(P0, P1, Q0)) intersect = true;
            if (std::abs(o2) < eps && on_segment(P0, P1, Q1)) intersect = true;
            if (std::abs(o3) < eps && on_segment(Q0, Q1, P0)) intersect = true;
            if (std::abs(o4) < eps && on_segment(Q0, Q1, P1)) intersect = true;
        }

        // ---------------------------------------
        // 5. Near-miss (endpoint very close) FIX
        // ---------------------------------------
        if (!intersect) {
            auto point_segment_dist2 = [](const Vec2& p, const Vec2& a, const Vec2& b) {
                Vec2 ab = b - a;
                double denom = dot(ab, ab);
                if (denom < 1e-12) return dot(p - a, p - a);

                double t = dot(p - a, ab) / denom;
                t = std::max(0.0, std::min(1.0, t));
                Vec2 proj = a + t * ab;
                return dot(p - proj, p - proj);
                };

            if (point_segment_dist2(P1, Q0, Q1) < eps_dist ||
                point_segment_dist2(Q0, P0, P1) < eps_dist ||
                point_segment_dist2(Q1, P0, P1) < eps_dist)
            {
                intersect = true;
            }
        }

        // ---------------------------------------
        // 6. Final decision
        // ---------------------------------------
        if (intersect && manifold.precond_flip_edge(flip_maybe)) {
            if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
                std::cout << "flipped!" << "\n";
            }
            return flip_maybe;
        }
    }

    return std::nullopt;
}

/// If an edge plane would be crossed, returns that edge to be flipped
std::optional<HalfEdgeID> find_crossed_edge(
    const Manifold& manifold,
    const VertexID id,
    const Point& starting_pos,
    const Point& end_pos,
    const ReexpandOpts& opts)
{
    for (auto he : manifold.incident_halfedges(id)) {
        auto walker = manifold.walker(he);
        auto flip_maybe = walker.next().halfedge();
        auto v1 = manifold.positions[walker.vertex()];
        auto v2 = manifold.positions[walker.next().vertex()];
        auto tangent = normalize(v1 - v2);
        if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
            std::cout << v1 << "\n";
            std::cout << v2 << "\n";
        }

        auto normal = CGLA::normalize(CGLA::cross((starting_pos - v1), tangent));

        // normal of the plane
        Vec3 binormal = CGLA::cross(tangent, normal);
        // point on the plane
        const Point& p0 = (v1 + v2) * 0.5;
        const double r = (v1 - v2).length() * 0.5;

        // starting point in our line
        const Point& l0 = starting_pos;
        const Vec3 l = normalize(end_pos - starting_pos);

        auto den = dot(l, binormal);
        if (std::abs(den) < 1e-8) {
            continue;
        } // parallel case
        auto len = (end_pos - starting_pos).length();
        auto d = dot((p0 - l0), binormal) / den;

        auto p = l0 + d * l;
        auto distance = (p - p0).length();
        if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
            std::cout << "d  : " << d << "\n";
            std::cout << "len: " << len << "\n";
        }
        if (distance > r) {
            continue;
        }

        if (d > 0 && len > (d * 0.90) && manifold.precond_flip_edge(flip_maybe)) {
            if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
                std::cout << "flipped!" << "\n";
            }
            return flip_maybe;
        }
    }
    return std::nullopt;
}

// Fast quick sort implementations from cppreference
// The reason why this is here is using std::sort results in memory corruption which probably has to do with
// incorrect NaN handling or some other poorly specified behavior in operator< for double.
namespace
{
    template <typename Comp = std::less<>>
    void quick_sort(std::forward_iterator auto first, std::forward_iterator auto last, Comp comp = std::less<>{})
    {
        if (first == last)
            return;

        auto pivot = *std::next(first, std::distance(first, last) / 2);
        auto middle1 = std::partition(first, last, [&](const auto& elem) {
            return comp(elem, pivot);
        });
        auto middle2 = std::partition(middle1, last, [&](const auto& elem) {
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

void fix_flipped_edge(Manifold& manifold, VertexID center_v, Vec3& normal) {
    bool isValid = false;

    while (!isValid) {
        bool isFlipped = false;
        for (auto h : manifold.incident_halfedges(center_v)) {
            if (manifold.walker(h).face() == InvalidFaceID && manifold.walker(h).opp().face() != InvalidFaceID) {
                if (!manifold.in_use(h))
                    continue;
                Point v1 = manifold.positions[manifold.walker(h).opp().next().vertex()];
                Point v_root = manifold.positions[manifold.walker(h).vertex()];
                Point v0 = manifold.positions[center_v];
                Vec3 cross_p = CGLA::cross(v0 - v_root, v1 - v_root);
                /*std::cout << normal << std::endl;
                std::cout << cross_p << std::endl;
                std::cout << "flip_case1" << std::endl;
                std::cout << center_v << std::endl;
                std::cout << (CGLA::dot(cross_p, normal) <= 0.) << std::endl;*/
                if (CGLA::dot(cross_p, normal) <= 0.) {
                    isFlipped = true;
                    manifold.remove_edge(h);
                    break;
                }
            }
            else if (manifold.walker(h).opp().face() == InvalidFaceID && manifold.walker(h).face() != InvalidFaceID) {
                if (!manifold.in_use(h))
                    continue;
                Point v1 = manifold.positions[manifold.walker(h).next().vertex()];
                Point v_root = manifold.positions[manifold.walker(h).vertex()];
                Point v0 = manifold.positions[center_v];
                Vec3 cross_p = CGLA::cross(v1 - v_root, v0 - v_root);
                /*std::cout << cross_p << std::endl;
                std::cout << "flip_case2" << std::endl;*/
                if (CGLA::dot(cross_p, normal) <= 0.) {
                    isFlipped = true;
                    manifold.remove_edge(h);
                    break;
                }
            }
        }
        if (!isFlipped)
            isValid = true;
    }
    return;
}

void showProgressBar(float progress) {
    int barWidth = 70;  // Width of the progress bar

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";  // \r returns to the beginning of the line
    std::cout.flush();  // Flush the output to show the progress
}

void reexpand_points(Manifold& manifold, const Collapse& collapse, const ReexpandOpts& opts)
{
    std::cout << "reexpanding" << std::endl;
    const auto& manifold_positions = manifold.positions;

    std::unordered_multimap<Point, VertexID, PointHash, PointEquals> point_to_manifold_ids;
    for (auto manifold_vid : manifold.vertices()) {
        auto pos = manifold_positions[manifold_vid];
        point_to_manifold_ids.emplace(pos, manifold_vid);
    }
    auto position_to_manifold_iter = [&](const Point& point) {
        auto [fst, snd] = point_to_manifold_ids.equal_range(point);
        return std::ranges::subrange(fst, snd) | std::views::values;
    };

    // insert latent point to stored latent position
    // update active point position to the stored coordinate

    double angle_threshold_cos = std::cos(opts.min_angle_threshold);
    // Now we need to consider two position candidates

    size_t expansion_failures = 0;
    size_t bad_expansions = 0;
    size_t flips = 0;
    int iteration = 0;
    int step = 0;
    int total_collapse = 0;
    for (const auto& collapse_iter : collapse.collapses | std::views::reverse) {
        total_collapse += collapse_iter.size();
    }
    std::cout << "Total reexpansion: " << total_collapse << std::endl;
    std::vector<HalfEdgeID> one_ring;
    std::vector<HalfEdgeID> circle;
    std::vector<HalfEdgeID> two_ring;

    /*int save_step1 = int(total_collapse * 0.2);
    int save_step2 = int(total_collapse * 0.5);
    int save_step3 = int(total_collapse * 0.8);*/


    for (const auto& collapse_iter : collapse.collapses | std::views::reverse) {
        iteration++;
        for (auto single_collapse : collapse_iter | std::views::reverse) {
            step++;
            /*if(step == save_step1)
                HMesh::obj_save("step_1.obj", manifold);
            if (step == save_step2)
                HMesh::obj_save("step_2.obj", manifold);
            if (step == save_step3)
                HMesh::obj_save("step_3.obj", manifold);*/

            //showProgressBar(float(step) / float(total_collapse));
            //std::cout << "Total reexpansion: " << total_collapse << std::endl;
            //std::cout << "current_step: " << step << std::endl;
            bool output_all = (step == 7);
            output_all = false;
            // find the manifold_ids for the active vertex

            const auto active_pos = single_collapse.active_point_coords;
            const auto latent_pos = single_collapse.latent_point_coords;
            const auto v_bar = single_collapse.v_bar;
            
            // FIXME: Debug info
            if (opts.debug_opts.debug_mask & RE_ITERATION) {
                std::cout << "--------------------------------\n";
                std::cout << "@iteration: " << iteration << "\n";
                std::cout << "active pos: " << active_pos << "\n";
                std::cout << "latent pos: " << latent_pos << "\n";
                std::cout << "combin pos: " << v_bar << "\n";
            }
            if (opts.debug_opts.debug_mask & RE_MARK_SPLITS) {
                manifold.add_face({active_pos, latent_pos, v_bar});
            }
            if (opts.debug_opts.debug_mask & RE_CROSSING_FLIP) {
                std::cout << "First flip:\n";
            }
            // repair local geometry maybe
            const auto manifold_ids = position_to_manifold_iter(v_bar);
            Vec3 v_bar_norm;
            bool isFailed = false;
            for (const auto id : manifold_ids) {
                v_bar_norm = manifold.normal(id);
                if (false) {
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
                // Fix crossed edge new
                else{
                    bool isBDRY_removed = false;
                    bool isValid = false;
                    int find_time = 0;

                    while (!isValid) {
                        isValid = true;
                        auto maybe1 = find_crossed_edge(manifold, id, active_pos, latent_pos, opts);
                        //auto maybe1 = find_crossed_edge_new(manifold, id, v_bar, active_pos, opts);
                        if (maybe1) {
                            // Decide remove or flip
                            if (manifold.walker(*maybe1).opp().face() == InvalidFaceID) {
                                // Intersecting boundary edge -- remove this boundary edge
                                manifold.remove_edge(*maybe1);
                                isBDRY_removed = true;
                            }
                            else {
                                if (manifold.precond_flip_edge(*maybe1))
                                    manifold.flip_edge(*maybe1);
                                else {
                                    //std::cout << "Can't flip" << std::endl;
                                    isFailed = true;
                                }
                                isValid = false;
                            }
                        }
                        auto maybe2 = find_crossed_edge(manifold, id, latent_pos, active_pos, opts);
                        //auto maybe2 = find_crossed_edge_new(manifold, id, v_bar, latent_pos, opts);
                        if (maybe2) {
                            // Decide remove or flip
                            if (manifold.walker(*maybe2).opp().face() == InvalidFaceID) {
                                // Intersecting boundary edge -- remove this boundary edge
                                manifold.remove_edge(*maybe2);
                                isBDRY_removed = true;
                            }
                            else {
                                if (manifold.precond_flip_edge(*maybe2)) {
                                    manifold.flip_edge(*maybe2);
                                }
                                else {
                                    //std::cout << "Can't flip" << std::endl;
                                    isFailed = true;
                                }
                                isValid = false;
                            }
                        }
                        find_time++;
                        if (find_time > 100) {
                            std::cout << "Failed" << std::endl;
                            isFailed = true;
                        }
                        if (isFailed)
                            break;
                    }

                    //std::cout << "crossed_edge fixed" << std::endl;

                    if (output_all)
                        HMesh::obj_save("C:/Users/ruicu/Desktop/SGP26/Results/debug/" +
                            std::to_string(step) + "_after_flip.obj", manifold);

                    manifold.positions[id] = active_pos;
                    if (output_all)
                        HMesh::obj_save("C:/Users/ruicu/Desktop/SGP26/Results/debug/" +
                            std::to_string(step) + "_after_reposition.obj", manifold);
                    if (false) {
                        if (v_bar_norm.length() < 1e-8 || std::isnan(v_bar_norm.length()) ||
                            std::isinf(v_bar_norm.length())) {
                            isFailed = true;
                            break;
                        }
                        fix_flipped_edge(manifold, id, v_bar_norm);
                    }

                    if (output_all)
                        HMesh::obj_save("C:/Users/ruicu/Desktop/SGP26/Results/debug/" +
                            std::to_string(step) + "_after_fix.obj", manifold);
                } 
            }
            
            if (isFailed) {
                expansion_failures++;
                continue;
            }
            // actually do the expansion
            //std::cout << "find edge pair" << std::endl;
            Split info;
            VertexID org_v;
            Vec3 org_norm;
            const auto new_vid = [&]() -> VertexID {
                // we want to get as close to 90 degrees as possible here
                for (const auto this_vert : manifold_ids) {
                    const auto candidate = find_edge_pair(manifold, this_vert, latent_pos, active_pos, opts,
                                                          angle_threshold_cos);
                    org_norm = manifold.normal(this_vert);
                    //const auto candidate = find_edge_pair(manifold, this_vert, latent_pos, active_pos);
                    
                    if (candidate.h_in != InvalidHalfEdgeID) {
                        const auto vnew = manifold.split_vertex(candidate.h_in, candidate.h_out);
                        GEL_ASSERT_NEQ(vnew, InvalidVertexID);
                        info = candidate;
                        org_v = this_vert;
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
            //std::cout << "find edge pair done" << std::endl;

            // Fix same side case: find the flipped face and newly added edge, remove them
            /*std::cout << info.isSameSide << std::endl;
            std::cout << manifold.walker(info.h_in).opp().vertex() << std::endl;
            std::cout << manifold.walker(info.h_out).vertex() << std::endl;
            */
            if (false) {
                bool isValid = false;
                while (!isValid) {
                    bool isRemove = false;
                    isValid = true;
                    for (const auto edge_id : manifold.incident_halfedges(org_v)) {
                        if (manifold.walker(edge_id).face() != InvalidFaceID) {
                            HalfEdgeID fan_edge = manifold.walker(edge_id).next().halfedge();
                            if (manifold.walker(fan_edge).opp().face() == InvalidFaceID) {
                                Vec3 vec1 = manifold.positions[manifold.walker(fan_edge).vertex()] -
                                    manifold.positions[manifold.walker(fan_edge).opp().vertex()];
                                Vec3 vec2 = manifold.positions[manifold.walker(fan_edge).next().vertex()] -
                                    manifold.positions[manifold.walker(fan_edge).vertex()];
                                if (CGLA::dot(CGLA::cross(vec1, vec2), org_norm) < 0.) {
                                    manifold.remove_edge(fan_edge);
                                    isRemove = true;
                                    break;
                                }
                                //std::cout << "Triangle " << manifold.walker(split_edge).next().vertex() << " removed" << std::endl;
                            }
                        }
                    }

                    for (const auto edge_id : manifold.incident_halfedges(new_vid)) {
                        if (manifold.walker(edge_id).face() != InvalidFaceID) {
                            HalfEdgeID fan_edge = manifold.walker(edge_id).next().halfedge();
                            if (manifold.walker(fan_edge).opp().face() == InvalidFaceID) {
                                Vec3 vec1 = manifold.positions[manifold.walker(fan_edge).vertex()] -
                                    manifold.positions[manifold.walker(fan_edge).opp().vertex()];
                                Vec3 vec2 = manifold.positions[manifold.walker(fan_edge).next().vertex()] -
                                    manifold.positions[manifold.walker(fan_edge).vertex()];
                                if (CGLA::dot(CGLA::cross(vec1, vec2), org_norm) < 0.) {
                                    manifold.remove_edge(fan_edge);
                                    isRemove = true;
                                    break;
                                }
                                //std::cout << "Triangle " << manifold.walker(split_edge).next().vertex() << " removed" << std::endl;
                            }
                        }
                    }
                    if (isRemove)
                        isValid = false;
                }
                
            }

            if (false) {
                if (info.isSameSide) {
                    // find the halfedge {v_old, v_new}
                    HalfEdgeID split_edge;
                    for (const auto edge_id : manifold.incident_halfedges(org_v)) {
                        if (manifold.walker(edge_id).vertex() == new_vid)
                            split_edge = edge_id;
                    }

                    // Check if both face align with original normal, if not, delete the new face and new edge
                    // First triangle
                    if (manifold.walker(split_edge).face() != InvalidFaceID) {
                        Vec3 v_split_edge = manifold.positions[manifold.walker(split_edge).vertex()] -
                            manifold.positions[manifold.walker(split_edge).opp().vertex()];
                        Vec3 v_next_edge = manifold.positions[manifold.walker(split_edge).next().vertex()] -
                            manifold.positions[manifold.walker(split_edge).vertex()];
                        if (CGLA::dot(CGLA::cross(v_next_edge, -v_split_edge), org_norm) < 0.) {
                            HalfEdgeID to_remove = manifold.walker(split_edge).next().halfedge();
                            HalfEdgeID to_remove_2 = manifold.walker(to_remove).next().halfedge();
                            if (manifold.walker(to_remove).opp().face() == InvalidFaceID)
                                manifold.remove_edge(to_remove);
                            else if (manifold.walker(to_remove_2).opp().face() == InvalidFaceID)
                                manifold.remove_edge(to_remove_2);
                        }
                        //std::cout << "Triangle " << manifold.walker(split_edge).next().vertex() << " removed" << std::endl;
                    }
                    // Second triangle
                    HalfEdgeID split_edge_opp = manifold.walker(split_edge).opp().halfedge();
                    if (manifold.walker(split_edge_opp).face() != InvalidFaceID) {
                        Vec3 v_split_edge = manifold.positions[manifold.walker(split_edge_opp).vertex()] -
                            manifold.positions[manifold.walker(split_edge_opp).opp().vertex()];
                        Vec3 v_next_edge = manifold.positions[manifold.walker(split_edge_opp).next().vertex()] -
                            manifold.positions[manifold.walker(split_edge_opp).vertex()];
                        if (CGLA::dot(CGLA::cross(v_next_edge, -v_split_edge), org_norm) < 0.) {
                            HalfEdgeID to_remove = manifold.walker(split_edge_opp).next().halfedge();
                            HalfEdgeID to_remove_2 = manifold.walker(to_remove).next().halfedge();
                            if (manifold.walker(to_remove).opp().face() == InvalidFaceID)
                                manifold.remove_edge(to_remove);
                            else if (manifold.walker(to_remove_2).opp().face() == InvalidFaceID)
                                manifold.remove_edge(to_remove_2);
                        }
                        //std::cout << "Triangle " << manifold.walker(split_edge).next().vertex() << " removed" << std::endl;

                    }
                }

                // Repair flipped triangles
                // Check org_v
                bool no_flip = false;
                int num_flip = 0;
                while (!no_flip) {
                    no_flip = true;
                    for (const auto h : manifold.incident_halfedges(org_v)) {
                        const auto& h_walker = manifold.walker(h);
                        if (h_walker.face() != InvalidFaceID) {
                            Point v1 = manifold.positions[h_walker.next().vertex()];
                            Point v_root = manifold.positions[h_walker.vertex()];
                            Point v0 = manifold.positions[org_v];
                            Vec3 cross_p = CGLA::cross(v1 - v_root, v0 - v_root);
                            if (CGLA::dot(cross_p, org_norm) <= 0. && manifold.precond_flip_edge(h)) {
                                manifold.flip_edge(h);
                                no_flip = false;
                                break;
                            }
                        }
                    }
                    num_flip++;
                    if (num_flip > 100) {
                        isFailed = true;
                        break;
                    }
                }
                if (isFailed) {
                    expansion_failures++;
                    continue;
                }
                
                // Check new_v
                no_flip = false;
                num_flip = 0;
                while (!no_flip) {
                    no_flip = true;
                    for (const auto h : manifold.incident_halfedges(new_vid)) {
                        const auto& h_walker = manifold.walker(h);
                        if (h_walker.face() != InvalidFaceID) {
                            Point v1 = manifold.positions[h_walker.next().vertex()];
                            Point v_root = manifold.positions[h_walker.vertex()];
                            Point v0 = manifold.positions[new_vid];
                            Vec3 cross_p = CGLA::cross(v1 - v_root, v0 - v_root);
                            if (CGLA::dot(cross_p, org_norm) <= 0. && manifold.precond_flip_edge(h)) {
                                manifold.flip_edge(h);
                                no_flip = false;
                                break;
                            }
                        }
                    }
                    num_flip++;
                    if (num_flip > 100) {
                        isFailed = true;
                        break;
                    }
                }

                if (isFailed) {
                    expansion_failures++;
                    continue;
                }

            }
            
            //std::cout << "post-processing done" << std::endl;

            // Update the point to manifold id map
            // we need to copy the data to avoid invalidating iterators when we mutate point_to_manifold_ids
            if (std::ranges::distance(manifold_ids) < 4) {
                InplaceVector<VertexID, 3> copy(manifold_ids.begin(), manifold_ids.end());
                for (auto id : copy) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
            } else {
                std::vector copy(manifold_ids.begin(), manifold_ids.end());
                for (auto id : copy) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
            }
            point_to_manifold_ids.emplace(latent_pos, new_vid);
            point_to_manifold_ids.erase(v_bar);

            //HMesh::obj_save("C:/Users/ruicu/Desktop/SGP26/Results/debug/" +
                //std::to_string(step) + "_after_exp.obj", manifold);
            //std::cout << new_vid << std::endl;
            // populate a bunch of vectors for the optimization phase
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

            //std::cout << "start optimizing" << std::endl;
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

            // perform optimization/refinement
            //std::cout << opts.refinement_iterations << std::endl;
            for (int i = 0; i < opts.refinement_iterations; ++i) {
                auto threshold = opts.refinement_angle_threshold;
                quick_sort(one_ring, cmp);
                for (HalfEdgeID h : one_ring | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, threshold);
                    //bool flipped = angle_flip_check_new(manifold, h);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
                quick_sort(circle, cmp);
                for (HalfEdgeID h : circle | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, threshold);
                    //bool flipped = angle_flip_check_new(manifold, h);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
                quick_sort(two_ring, cmp);
                for (HalfEdgeID h : two_ring | std::views::reverse) {
                    bool flipped = angle_flip_check(manifold, h, threshold);
                    //bool flipped = angle_flip_check_new(manifold, h);
                    if (flipped) {
                        manifold.flip_edge(h);
                    }
                    flips += flipped;
                }
            }

            //HMesh::obj_save("C:/Users/ruicu/Desktop/SGP26/Results/debug/" +
                //std::to_string(step) + ".obj", manifold);
            if (opts.debug_opts.debug_mask & RE_ERRORS) {
                const auto v_new_max_angle = one_ring_max_dihedral_angle(manifold, new_vid);
                const auto v_old_max_angle = one_ring_max_dihedral_angle(manifold, manifold_ids.front());
                if (v_new_max_angle > opts.debug_opts.angle_bad_threshold || v_old_max_angle > opts.debug_opts.
                    angle_bad_threshold) {
                    bad_expansions++;
                    if (opts.debug_opts.debug_mask & RE_ERRORS) {
                        auto normal = manifold.normal(new_vid);
                        manifold.add_face({
                            active_pos, latent_pos, v_bar + (active_pos - latent_pos).length() * normal
                        });
                        std::cout << "v_new max angle: " << v_new_max_angle << "\n";
                        std::cout << "v_old max angle: " << v_old_max_angle << "\n";
                        std::cout << "failed at iteration " << iteration << "\n";
                    }
                    if (opts.debug_opts.early_stop_at_error && opts.debug_opts.stop_at_error == 0) {
                        goto EXIT;
                    }
                    if (opts.debug_opts.stop_at_error > 0 && opts.debug_opts.stop_at_error == bad_expansions) {
                        goto EXIT;
                    }
                }
            }

            if (opts.debug_opts.stop_at_iteration > 0 && iteration == opts.debug_opts.stop_at_iteration) {
                std::cout << "stopped early at " << iteration << "\n";
                goto EXIT;
            }
        }
        // DEBUG
        //HMesh::obj_save("iteration_" + std::to_string(iteration) + ".obj", manifold);
    }
EXIT:
    std::cout << "flips: " << flips << "\n";
    std::cerr << "failures: " << expansion_failures << "\n";
}
}
