//
// Created by Cem Akarsubasi on 9/10/25.
//

#include "Collapse.h"

#include <unordered_map>
#include <GEL/Util/RangeTools.h>

namespace HMesh::RSR
{
struct Split {
    HMesh::HalfEdgeID h_in;
    HMesh::HalfEdgeID h_out;
    std::optional<HalfEdgeID> flip = std::nullopt;
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

Split find_edge_pair(const HMesh::Manifold& m, const HMesh::VertexID center_idx, const Vec3& v_new_position,
    const Vec3& v_old_position, const ReexpandOptions& opts) {
    // this is getting too complicated
    const auto print_hedge = [&m](HalfEdgeID he) {
        auto v_from = m.positions[m.walker(he).opp().vertex()];
        auto v_to = m.positions[m.walker(he).vertex()];
        std::cout << v_from << " -> " << v_to << " (" << (v_to - v_from) << ")\n";
    };

    const auto triangle_from_half_edge_orig = [&m](const Point& origin, const Walker& w) -> Triangle {
        auto v1 = w.vertex();
        auto v2 = w.next().vertex();
        return Triangle{origin, m.positions[v1], m.positions[v2]};
    };

    const auto triangle_from_half_edge = [&m](const Walker& w) -> Triangle {
        auto v1 = w.vertex();
        auto v2 = w.next().vertex();
        auto v3 = w.next().next().vertex();
        return Triangle{m.positions[v1], m.positions[v2], m.positions[v3]};
    };

    // Optimize the dihedral value
    const auto dihedral_one_ring = [](const Triangle& t1, const Triangle& t2) -> double {
        auto d = optimize_dihedral(t1.normal(), t2.normal());
        return std::pow(d, 1.0) / std::min(t1.area(), t2.area()); // FIXME: figure this out
    };

    const auto dihedral_two_ring = [&m](const Triangle& t, const FaceID f) -> double {
        if (f == InvalidFaceID) return 0.0;
        const auto d = optimize_dihedral(t.normal(), m.normal(f));
        return std::pow(d, 1.0) / std::min(m.area(f), t.area());
    };

    // find the edge flip candidate now
    struct FlipCandidate {
        HalfEdgeID candidate; // the actual candidate to be flipped
        HalfEdgeID he1; // from the previous two
        HalfEdgeID he2; // from the previous two
        Point opposing_vertex; // required to calculate dihedrals
    };
    auto find_flip_candidate = [&m](HalfEdgeID he1, HalfEdgeID he2) -> std::optional<FlipCandidate> {
        auto v_target = m.walker(he2).vertex();
        auto walker = m.walker(he1);
        auto v_orig = walker.opp().vertex();
        for (; !walker.full_circle(); walker = walker.circulate_vertex_ccw()) {
            if (walker.vertex() == v_target) {
                auto candidate = walker.halfedge();
                Point opposing_vertex;
                if (walker.next().vertex() == v_orig) {
                    opposing_vertex = m.positions[walker.opp().next().vertex()];
                } else {
                    opposing_vertex = m.positions[walker.next().vertex()];
                    std::swap(he1, he2);
                }
                return FlipCandidate{
                    .candidate = candidate,
                    .he1 = he1,
                    .he2 = he2,
                    .opposing_vertex = opposing_vertex
                };
            }
        }
        return std::nullopt;
    };

    struct ExpandResult {
        double score = INFINITY;
        double max_angle = INFINITY;
        std::vector<Triangle> triangles; // FIXME: debug time
        // TODO: maybe add an optional "flip"
        std::optional<HalfEdgeID> flip = std::nullopt;
    };
    const auto expand_score = [&](
        HalfEdgeID h_in_opp,
        HalfEdgeID h_out,
        const Point& v_new_position,
        const Point& v_old_position,
        const std::optional<FlipCandidate>& flip_candidate) -> ExpandResult {
        const auto walker_out = m.walker(h_out);
        const auto walker_in_opp = m.walker(h_in_opp);
        GEL_ASSERT(flip_candidate.has_value());
        const auto v_h_out = m.positions[walker_out.vertex()];
        const auto v_h_in = m.positions[walker_in_opp.vertex()];

        // instead of doing this unscalable stupidity, let's try to be smart and perform a rotation through all of the
        // affected triangles. We basically perform a sliding window and construct the right triangle by
        // either passing v_new or v_old as the third point. The order of the window also matters a lot

        const auto tri_center_in = Triangle(v_old_position, v_new_position, v_h_in);
        const auto tri_center_out = Triangle(v_old_position, v_h_out, v_new_position);

        double two_ring = 0.0;
        auto one_ring_iterator = [&]() {
            std::vector<Triangle> triangles;

            auto walker_prev = walker_out;
            auto walker_next = walker_prev.circulate_vertex_ccw();
            // from out towards in counterclockwise
            triangles.push_back(tri_center_out);
            while (walker_prev.halfedge() != walker_in_opp.halfedge()) {
                auto p2 = m.positions[walker_prev.vertex()];
                auto p3 = m.positions[walker_next.vertex()];
                Triangle tri = {v_new_position, p2, p3};

                if (walker_prev.halfedge() == flip_candidate->he1 && walker_next.halfedge() == flip_candidate->he2) {
                    FaceID opposing_face = walker_prev.next().opp().face();
                    two_ring = dihedral_two_ring(tri, opposing_face);
                }

                // to consider the two ring dihedrals, we need to get the triangles from the opposite edges.
                // none of the triangles are affected by the expansion, so we can just fetch them from the manifold directly
                //


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

                Triangle tri = {v_old_position, p2, p3};

                triangles.push_back(tri);
                walker_prev = walker_prev.circulate_vertex_ccw();
                walker_next = walker_next.circulate_vertex_ccw();
            }
            return triangles;
        };
        // TODO
        auto one_ring_iterator_flipped_edge = [&]() {
            std::vector<Triangle> triangles;

            auto walker_prev = walker_out;
            auto walker_next = walker_prev.circulate_vertex_ccw();
            // from out towards in counterclockwise
            triangles.push_back(tri_center_out);
            while (walker_prev.halfedge() != walker_in_opp.halfedge()) {
                if (walker_prev.halfedge() == flip_candidate->he1 && walker_next.halfedge() == flip_candidate->he2) {
                    // TODO
                    auto p2 = m.positions[walker_prev.vertex()];
                    auto p3 = m.positions[walker_next.vertex()];
                    Triangle tri1 = {v_new_position, flip_candidate->opposing_vertex, p2};
                    Triangle tri2 = {v_new_position, p3, flip_candidate->opposing_vertex};
                    triangles.push_back(tri1);
                    triangles.push_back(tri2);
                } else {
                    auto p2 = m.positions[walker_prev.vertex()];
                    auto p3 = m.positions[walker_next.vertex()];

                    Triangle tri = {v_new_position, p2, p3};

                    // to consider the two ring dihedrals, we need to get the triangles from the opposite edges.
                    // none of the triangles are affected by the expansion, so we can just fetch them from the manifold directly
                    //FaceID opposing_face = walker_prev.next().opp().face();
                    //two_ring += dihedral_two_ring(tri, opposing_face);

                    triangles.push_back(tri);
                }
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

                Triangle tri = {v_old_position, p2, p3};

                triangles.push_back(tri);
                walker_prev = walker_prev.circulate_vertex_ccw();
                walker_next = walker_next.circulate_vertex_ccw();
            }
            return triangles;
        };

        // h_in_opp and h_out are unique, making this sound
        // TODO: get rid of this
        //double two_ring = 0.0;
        auto triangles = one_ring_iterator();

        auto dihedral0 = dihedral_one_ring(tri_center_in, tri_center_out);

        auto calculate_dihedrals = [&](auto& iter) {
            double total_dihedral = 0.0;
            double max_angle = 0.0;
            auto shifted = Util::Ranges::shifted_wrapping(iter, 1);
            for (auto [tri1, tri2] : Util::Ranges::zip(iter, shifted)) {
                const auto one_ring = dihedral_one_ring(tri1, tri2);
                max_angle = std::max(max_angle, one_ring);
                total_dihedral += one_ring;
            }
            max_angle = std::max(max_angle, dihedral0);
            return std::make_pair(total_dihedral + dihedral0, max_angle);
        };
        auto [total_dihedral, max_angle] = calculate_dihedrals(triangles);
        auto triangles2 = one_ring_iterator_flipped_edge();
        auto [total_dihedral_s, max_angle_s] = calculate_dihedrals(triangles2);

        auto valency_cost = 0; //optimize_valency(m, h_out, h_in_opp, opts);
        auto total_cost = total_dihedral + dihedral0 + valency_cost;
        //total_dihedral_s = 100000.0;
        if (total_dihedral < total_dihedral_s) {
            return ExpandResult{
                .score = total_dihedral + two_ring,
                .max_angle = max_angle,
                .triangles = std::move(triangles),
            };
        } else {
            return ExpandResult{
                .score = total_dihedral_s,
                .max_angle = max_angle_s,
                .triangles = std::move(triangles2),
                .flip = flip_candidate->candidate
            };
        }
    };

    // let's do it the dumb way for once
    std::vector<HalfEdgeID> half_edges;
    HMesh::circulate_vertex_ccw(m, center_idx, [&](HalfEdgeID he) {
        half_edges.emplace_back(he);
    });


    struct HalfEdgeHelper {
        HalfEdgeID he;
        double prod = INFINITY;
    };
    struct HalfEdgeHelperComp {
        bool operator()(const HalfEdgeHelper& lhs, const HalfEdgeHelper& rhs) const
        {
            return lhs.prod < rhs.prod;
        }
    };

    double min_score = INFINITY;
    double max_angle = INFINITY;
    HalfEdgeID h_in_opp;
    HalfEdgeID h_out;
    //std::optional<HalfEdgeID> flip_candidate;
    std::vector<Triangle> triangles;

    // TODO: figure out the flip candidate
    Util::InplaceVector<HalfEdgeHelper, 3> hedges;
    //std::vector<HalfEdgeHelper> hedges;
    auto edge = CGLA::normalize(v_new_position - v_old_position);
    for (auto h : half_edges) {
        auto v = m.positions[m.walker(h).vertex()];
        auto e = CGLA::normalize(v - v_new_position);
        auto prod = std::abs(CGLA::dot(e, edge) - 1.0);
        hedges.emplace_back(h, prod);
        // TODO: bubble sort would be good here
        std::ranges::sort(hedges, HalfEdgeHelperComp());
        if (hedges.size() == 3) {
            hedges.pop_back();
        }
    }
    GEL_ASSERT_EQ(hedges.size(), 2);
    auto flip_candidate = find_flip_candidate(hedges[0].he, hedges[1].he);
    std::optional<HalfEdgeID> flip = std::nullopt;
    for (auto h1 : half_edges) {
        for (auto h2 : half_edges) {
            if (h1 == h2) {
                continue;
            }
            auto score = expand_score(h1, h2, v_new_position, v_old_position, flip_candidate);
            if (score.score < min_score) {
                min_score = score.score;
                max_angle = score.max_angle;
                triangles = std::move(score.triangles);
                flip = score.flip;
                h_in_opp = h1;
                h_out = h2;
            }
        }
    }
    if (flip.has_value()) {
        for (auto tri: triangles) {
            std::cout << tri << "\n";
        }
        std::cout << "------------";
    }


    return Split{m.walker(h_in_opp).opp().halfedge(), h_out, flip};
}

auto collapse_points(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
    const CollapseOpts& options) -> Collapse {
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

double optimize_dihedral(const Vec3& n1, const Vec3& n2) {
    const auto angle = CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0;
    return std::abs(angle);
    //const auto angle_cos = std::abs(CGLA::dot(n1, n2)) / (CGLA::length(n1) * CGLA::length(n2));
    //GEL_ASSERT_FALSE(std::isnan(angle_cos));
    //const auto angle = std::acos(angle_cos);
    //GEL_ASSERT_FALSE(std::isnan(angle));
    //return angle;
    //return std::abs(angle - 1.0);
}

double optimize_angle(const Vec3& p1, const Vec3& p2, const Vec3& p3, const ReexpandOptions& opts) {
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

double optimize_valency(const HMesh::Manifold& m, const HMesh::HalfEdgeID h_out, const HMesh::HalfEdgeID h_in_opp,
    const ReexpandOptions& opts) {
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

auto reexpand_points(HMesh::Manifold& manifold, Collapse2&& collapse, const ReexpandOptions& opts) -> void {
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
                if (candidate.flip.has_value()) {
                    manifold.flip_edge(*candidate.flip);
                }
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
            } else {
                failures++;
            }
        }
    }
    std::cerr << "failures: " << failures << "\n";
}

auto decimate(const Manifold& manifold, double factor) -> Manifold {
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

auto decimate_reexpand(const Manifold& manifold, double factor) -> Manifold {
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

    reexpand_points(m, std::move(collapse));

    return m;
}
}
