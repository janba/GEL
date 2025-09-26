//
// Created by Cem Akarsubasi on 9/10/25.
//

#include "Collapse.h"

#include <unordered_map>
#include <GEL/Util/RangeTools.h>

namespace HMesh::RSR
{
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
    const auto angle = CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0;
    return std::abs(angle);
}

// returns 0 for an equilateral triangle
double optimize_angle(const Vec3& p1, const Vec3& p2, const Vec3& p3, const ReexpandOptions& opts)
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
double optimize_valency(const HMesh::Manifold& m, const HMesh::HalfEdgeID h_out, const HMesh::HalfEdgeID h_in_opp,
                        const ReexpandOptions& opts)
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
    double shared_length;
    Vec3 normal;

    Triangle(const Vec3& p1, const Vec3& p2, const Vec3& p3, double shared_length) :
        shared_length(shared_length), normal(triangle_normal(p1, p2, p3))
    {}
};

Split find_edge_pair(const Manifold& m, const VertexID center_idx, const Vec3& v_new_position,
                     const Vec3& v_old_position, const ReexpandOptions& opts)
{
    // this is getting too complicated
    const auto print_hedge = [&m](HalfEdgeID he) {
        auto v_from = m.positions[m.walker(he).opp().vertex()];
        auto v_to = m.positions[m.walker(he).vertex()];
        std::cout << v_from << " -> " << v_to << " (" << (v_to - v_from) << ")\n";
    };

    // Optimize the dihedral value
    const auto dihedral_one_ring = [](const Triangle& t1, const Triangle& t2) -> double {
        auto d = optimize_dihedral(t1.normal, t2.normal);
        auto edge_length = t1.shared_length; // shared_edge_length(t1, t2);
        return d * edge_length;
    };

    std::vector<Triangle> triangles_buffer;

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
        const auto v_h_out = m.positions[walker_out.vertex()];
        const auto v_h_in = m.positions[walker_in_opp.vertex()];

        // instead of doing this unscalable stupidity, let's try to be smart and perform a rotation through all of the
        // affected triangles. We basically perform a sliding window and construct the right triangle by
        // either passing v_new or v_old as the third point. The order of the window also matters a lot
        const double in_len = (v_old_position - v_h_in).length();
        const double out_len = (v_new_position - v_h_out).length();
        const auto tri_center_in = Triangle(v_old_position, v_new_position, v_h_in, in_len);
        const auto tri_center_out = Triangle(v_new_position, v_old_position, v_h_out, out_len);

        auto one_ring_iterator = [&]() {
            triangles_buffer.clear();
            //std::vector<Triangle> triangles;

            auto walker_prev = walker_out;
            auto walker_next = walker_prev.circulate_vertex_ccw();
            // from out towards in counterclockwise
            triangles_buffer.push_back(tri_center_out);
            while (walker_prev.halfedge() != walker_in_opp.halfedge()) {
                const auto& p2 = m.pos(walker_prev.vertex());
                const auto& p3 = m.pos(walker_next.vertex());
                const auto shared_length = (v_new_position - p3).length();
                Triangle tri = {v_new_position, p2, p3, shared_length};
                // to consider the two ring dihedrals, we need to get the triangles from the opposite edges.
                // none of the triangles are affected by the expansion, so we can just fetch them from the manifold directly

                triangles_buffer.push_back(tri);
                walker_prev = walker_prev.circulate_vertex_ccw();
                walker_next = walker_next.circulate_vertex_ccw();
            }
            triangles_buffer.push_back(tri_center_in);
            // from in towards out counterclockwise
            walker_prev = walker_in_opp;
            walker_next = walker_in_opp.circulate_vertex_ccw();
            while (walker_prev.halfedge() != walker_out.halfedge()) {
                const auto& p2 = m.pos(walker_prev.vertex());
                const auto& p3 = m.pos(walker_next.vertex());
                const auto shared_length = (v_old_position - p3).length();
                Triangle tri = {v_old_position, p2, p3, shared_length};

                triangles_buffer.push_back(tri);
                walker_prev = walker_prev.circulate_vertex_ccw();
                walker_next = walker_next.circulate_vertex_ccw();
            }
            //return triangles;
        };

        // h_in_opp and h_out are unique, making this sound
        //auto triangles = one_ring_iterator();
        one_ring_iterator();

        const auto dihedral0 = dihedral_one_ring(tri_center_in, tri_center_out);
        auto calculate_dihedrals = [&](auto& iter) {
            double total_dihedral = 0.0;
            double max_angle = 0.0;
            // Not any faster, we need to do less work or cut out
            for (auto i = 0; i < triangles_buffer.size(); i++) {
                auto& tri1 = triangles_buffer[i];
                auto& tri2 = triangles_buffer[(i+1)%triangles_buffer.size()];
                const auto one_ring = dihedral_one_ring(tri1, tri2);
                max_angle = std::max(max_angle, one_ring);
                total_dihedral += one_ring;
            }
            //auto length = std::ranges::distance(iter);
            //auto repeated = Util::Ranges::repeat_range(iter | std::views::all) | std::views::take(length + 1);
            //auto shifted = Util::Ranges::shifted_wrapping(repeated, 1);
            //for (auto [tri1, tri2] : Util::Ranges::zip(repeated, shifted)) {
            //    static_assert(std::same_as<decltype(tri1), Triangle&>);
            //    const auto one_ring = dihedral_one_ring(tri1, tri2);
            //    max_angle = std::max(max_angle, one_ring);
            //    total_dihedral += one_ring;
            //}
            max_angle = std::max(max_angle, dihedral0);
            return std::make_pair(total_dihedral + dihedral0, max_angle);
        };
        auto [total_dihedral, max_angle] = calculate_dihedrals(triangles_buffer);
        return ExpandResult{
            .score = total_dihedral,
            .max_angle = max_angle,
        };
    };

    double min_score = INFINITY;
    double max_angle = INFINITY;
    HalfEdgeID h_in_opp;
    HalfEdgeID h_out;

    for (auto h1 : m.incident_halfedges(center_idx)) {
        for (auto h2 : m.incident_halfedges(center_idx)) {
            if (h1 == h2) {
                continue;
            }
            auto score = expand_score(h1, h2, v_new_position, v_old_position);
            if (score.score < min_score) {
                min_score = score.score;
                max_angle = score.max_angle;
                h_in_opp = h1;
                h_out = h2;
            }
        }
    }

    return Split{m.walker(h_in_opp).opp().halfedge(), h_out};
}

auto dihedral_from_half_edge(const Manifold& m, const HalfEdgeID h) -> double
{
    auto walker = m.walker(h);
    auto sp1 = m.positions[walker.vertex()];
    auto t1p3 = m.positions[walker.next().vertex()];
    auto sp2 = m.positions[walker.opp().vertex()];
    auto t2p3 = m.positions[walker.opp().next().vertex()];

    auto n1 = triangle_normal(sp1, t1p3, sp2);
    auto n2 = triangle_normal(sp2, t2p3, sp1);
    return m.length(h) * optimize_dihedral(n1, n2);
};

auto four_dihedrals(const Manifold& manifold, HalfEdgeID he) -> double
{
    const auto walker = manifold.walker(he);
    const auto he_n = walker.next().halfedge();
    const auto he_p = walker.prev().halfedge();
    const auto he_opp_n = walker.opp().next().halfedge();
    const auto he_opp_p = walker.opp().prev().halfedge();
    return dihedral_from_half_edge(manifold, he) +
        dihedral_from_half_edge(manifold, he_n) +
        dihedral_from_half_edge(manifold, he_p) +
        dihedral_from_half_edge(manifold, he_opp_n) +
        dihedral_from_half_edge(manifold, he_opp_p);
}

auto maybe_flip(Manifold& manifold, HalfEdgeID he, double dihedral_threshold = 1.5) -> bool
{
    if (!manifold.precond_flip_edge(he))
        return false;

    constexpr auto angle_threshold = 0.25;
    //constexpr auto dihedral_threshold = 1.5;
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
    const auto current_dihedral_cost = four_dihedrals(manifold, he);

    manifold.flip_edge(he);
    const auto new_dihedral_cost = four_dihedrals(manifold, he);
    if (new_dihedral_cost > current_dihedral_cost * dihedral_threshold) {
        manifold.flip_edge(he);
        return false;
    }

    if (angle1 < 0.0 || angle2 < 0.0 || angle3 < 0.0 || angle4 < 0.0) {
        manifold.flip_edge(he);
        return false;
        //return false;
    }
    return true;
    //if (((angle1 + angle2 < threshold) || (angle3 + angle4 < threshold)) && manifold.precond_flip_edge(he))
    //    manifold.flip_edge(he);
    if (((angle1 < angle_threshold && angle2 < angle_threshold) || (angle3 < angle_threshold && angle4 <
        angle_threshold))) {
        //manifold.flip_edge(he);
        return true;
    }
    manifold.flip_edge(he);
    return false;
}

auto collapse_points(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
                     const CollapseOpts& options) -> Collapse
{
    // set of connectivity information
    auto indices = [&vertices] {
        std::vector<NodeID> temp(vertices.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    Util::ImmediatePool pool;
    QuadraticCollapseGraph graph;

    // initialize graph
    for (auto i = 0UL; i < vertices.size(); ++i) {
        graph.add_node(vertices[i], normals[i]);
    }


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


    for (size_t iter = 0; iter < options.max_iterations; ++iter) {
        // TODO: stricter checking
        const size_t max_collapses = vertices.size() * std::pow(0.5, iter) * options.reduction_per_iteration;
        Collapse::ActivityMap activity_map;

        size_t count = 0;
        while (count < max_collapses) {
            count++;

            auto [active, latent, active_point_coords, latent_point_coords, v_bar] = graph.collapse_one();

            activity_map.insert(active, latent, active_point_coords, latent_point_coords, v_bar);
        }
        collapse.insert_collapse(activity_map);
    }
    //std::cout << "collapsed: " << count << std::endl;
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

auto reexpand_points(Manifold& manifold, Collapse2&& collapse, const ReexpandOptions& opts) -> void
{
    std::cout << "reexpanding\n";
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
    size_t flips = 0;
    for (const auto& collapse_iter : collapse.collapses | std::views::reverse) {
        for (auto single_collapse : collapse_iter | std::views::reverse) {
            // find the manifold_ids for the active vertex
            const auto active_pos = single_collapse.active_point_coords;
            const auto latent_pos = single_collapse.latent_point_coords;
            const auto v_bar = single_collapse.v_bar;
            const auto manifold_ids = position_to_manifold_iter(v_bar);
            for (const auto id : manifold_ids) {
                manifold.positions[id] = active_pos;
            }
            if (const auto new_vid = try_two_edge_expand(manifold_ids, latent_pos, active_pos); new_vid !=
                InvalidVertexID) {
                manifold.positions[new_vid] = latent_pos;
                for (const auto id : manifold_ids) {
                    point_to_manifold_ids.emplace(active_pos, id);
                }
                point_to_manifold_ids.emplace(latent_pos, new_vid);
                point_to_manifold_ids.erase(v_bar);

                // TODO: Look for edge flips
                // Our primary problem was on the actual ring around the

                // for (edge: one_ring)
                //   if ( dihedral angle is low AND minimum angle gain is positive AND flip is valid)
                //       perform a flip
                std::vector<HalfEdgeID> one_ring;
                std::vector<HalfEdgeID> circle;
                circulate_vertex_ccw(manifold, new_vid, [&](Walker& w) {
                    HalfEdgeID one_ring_he = w.halfedge();
                    one_ring.push_back(one_ring_he);
                    HalfEdgeID circle_he = w.next().halfedge();
                    circle.push_back(circle_he);
                });
                for (HalfEdgeID h : one_ring) {
                    flips += maybe_flip(manifold, h, 2.5);
                }
                for (HalfEdgeID h : circle) {
                    flips += maybe_flip(manifold, h, 2.5);
                }
                for (HalfEdgeID h : one_ring) {
                    flips += maybe_flip(manifold, h, 1.5);
                }
                for (HalfEdgeID h : circle) {
                    flips += maybe_flip(manifold, h, 1.5);
                }
                for (HalfEdgeID h : one_ring) {
                    flips += maybe_flip(manifold, h, 1.0);
                }
                for (HalfEdgeID h : circle) {
                    flips += maybe_flip(manifold, h, 1.0);
                }
            } else {
                failures++;
            }
        }
    }
    std::cout << "flips: " << flips << "\n";
    std::cerr << "failures: " << failures << "\n";
}

struct ExpansionChoice {
    /// The capacity on this one should determine how many edge pairs we even check
    Util::InplaceVector<HalfEdgeID, 3> to_flip;
};


auto decimate(const Manifold& manifold, double factor) -> Manifold
{
    if (factor >= 1.0 || factor < 0.0) {
        throw std::runtime_error("Invalid factor");
    }

    QuadraticCollapseGraph graph;
    // insert QEM for every point
    for (const auto vertex_id : manifold.vertices()) {
        graph.add_node(manifold.positions[vertex_id], manifold.normal(vertex_id));
    }
    // insert every edge into a queue
    for (const auto edge_id : manifold.halfedges()) {
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
    for (NodeID id : graph.node_ids()) {
        auto p1 = graph.m_vertices.at(id).position;
        vertices.push_back(p1);
        for (NodeID neighbor : graph.neighbors_lazy(id)) {
            if (id < neighbor) {
                for (NodeID third : graph.shared_neighbors(id, neighbor)) {
                    if (neighbor < third) {
                        auto n1 = graph.m_vertices.at(id).normal;
                        auto n2 = graph.m_vertices.at(neighbor).normal;
                        auto n3 = graph.m_vertices.at(third).normal;

                        auto p2 = graph.m_vertices.at(neighbor).position;
                        auto p3 = graph.m_vertices.at(third).position;
                        auto e1 = p2 - p1;
                        auto e2 = p3 - p1;
                        auto n = CGLA::cross(e1, e2);
                        if (CGLA::dot((n1 + n2 + n3) * (1. / 3.), n) > 0.0) {
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

auto decimate_reexpand(const Manifold& manifold, double factor) -> Manifold
{
    if (factor >= 1.0 || factor < 0.0) {
        throw std::runtime_error("Invalid factor");
    }

    QuadraticCollapseGraph graph;
    // insert QEM for every point
    for (const auto vertex_id : manifold.vertices()) {
        graph.add_node(manifold.positions[vertex_id], manifold.normal(vertex_id));
    }
    // insert every edge into a queue
    for (const auto edge_id : manifold.halfedges()) {
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
        collapse.collapses[0].emplace_back(single_collapse.active_point_coords, single_collapse.latent_point_coords,
                                           single_collapse.v_bar);
    }

    // create a new manifold from the collapsed graph
    std::vector<Vec3> vertices;
    vertices.reserve(graph.no_nodes());
    std::vector<NodeID> indices;
    indices.reserve(graph.no_edges() * 2);
    for (NodeID id : graph.node_ids()) {
        auto p1 = graph.m_vertices.at(id).position;
        vertices.push_back(p1);
        for (NodeID neighbor : graph.neighbors_lazy(id)) {
            if (id < neighbor) {
                for (NodeID third : graph.shared_neighbors(id, neighbor)) {
                    if (neighbor < third) {
                        auto n1 = graph.m_vertices.at(id).normal;
                        auto n2 = graph.m_vertices.at(neighbor).normal;
                        auto n3 = graph.m_vertices.at(third).normal;

                        auto p2 = graph.m_vertices.at(neighbor).position;
                        auto p3 = graph.m_vertices.at(third).position;
                        auto e1 = p2 - p1;
                        auto e2 = p3 - p1;
                        auto n = CGLA::cross(e1, e2);
                        if (CGLA::dot((n1 + n2 + n3) * (1. / 3.), n) > 0.0) {
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

    reexpand_points(m, std::move(collapse), ReexpandOptions());

    return m;
}
}
