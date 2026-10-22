/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/LRsR.h>

#include <GEL/Geometry/Delaunay3D.h>
#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/cleanup.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace HMesh::LRsR {
namespace {

using Vec3 = CGLA::Vec3d;

struct TriKey {
    size_t a = 0;
    size_t b = 0;
    size_t c = 0;
    bool operator==(const TriKey& other) const
    {
        return a == other.a && b == other.b && c == other.c;
    }
};

struct TriKeyHash {
    size_t operator()(const TriKey& tri) const noexcept
    {
        size_t h = tri.a;
        h ^= tri.b + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= tri.c + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

TriKey tri_key(size_t i, size_t j, size_t k)
{
    if (i > j) std::swap(i, j);
    if (j > k) std::swap(j, k);
    if (i > j) std::swap(i, j);
    return {i, j, k};
}

double median_neighbor_spacing(const std::vector<Vec3>& points)
{
    if (points.size() < 2)
        return 1.0;
    Geometry::Tree tree;
    for (size_t i = 0; i < points.size(); ++i) {
        const Vec3& point = points[i];
        if (std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]))
            tree.insert(point, i);
    }
    tree.build();

    std::vector<double> nearest;
    nearest.reserve(points.size());
    std::vector<Geometry::Tree::IdRecord> hits;
    for (size_t i = 0; i < points.size(); ++i) {
        hits.clear();
        tree.m_closest_ids(2, points[i], std::numeric_limits<double>::infinity(), hits);
        double best = 0.0;
        for (const Geometry::Tree::IdRecord& hit : hits) {
            if (hit.v == i || !(hit.d > 1e-24) || !std::isfinite(hit.d))
                continue;
            const double distance = std::sqrt(hit.d);
            if (best == 0.0 || distance < best)
                best = distance;
        }
        if (best > 0.0)
            nearest.push_back(best);
    }
    if (nearest.empty())
        return 1.0;
    const size_t mid = nearest.size() / 2;
    std::nth_element(nearest.begin(), nearest.begin() + static_cast<std::ptrdiff_t>(mid), nearest.end());
    if (nearest.size() % 2 == 1)
        return nearest[mid];
    const double upper = nearest[mid];
    const double lower = *std::max_element(nearest.begin(), nearest.begin() + static_cast<std::ptrdiff_t>(mid));
    return 0.5 * (lower + upper);
}

double triangle_area(const Vec3& a, const Vec3& b, const Vec3& c)
{
    return 0.5 * CGLA::length(CGLA::cross(b - a, c - a));
}

// Perimeter times area, divided by the area projected along `normal`.
// The projected area is the shadow on the plane perpendicular to that normal.
double projected_perimeter_score(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& normal)
{
    const Vec3 cross = CGLA::cross(b - a, c - a);
    const double cross_length = cross.length();
    const double normal_length = normal.length();
    if (!(cross_length > 0.0) || !(normal_length > 0.0)
        || !std::isfinite(cross_length) || !std::isfinite(normal_length))
        return std::numeric_limits<double>::infinity();
    const double projected_area = 0.5 * std::abs(CGLA::dot(cross, normal)) / normal_length;
    if (!(projected_area > 0.0) || !std::isfinite(projected_area))
        return std::numeric_limits<double>::infinity();
    const double area = 0.5 * cross_length;
    const double perimeter = CGLA::length(b - a) + CGLA::length(c - b) + CGLA::length(a - c);
    if (!std::isfinite(perimeter))
        return std::numeric_limits<double>::infinity();
    return perimeter * perimeter;// * area / projected_area;
}

// Both endpoints fit in 32 bits. reconstruct rejects a larger cloud.
uint64_t undirected_id(size_t a, size_t b)
{
    if (a > b)
        std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
}

// A triangle this close to the plane of the edge and its normal is neither side.
constexpr double kSideEpsilon = 1e-8;

bool finite_vector(const Vec3& value)
{
    return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

// Triangles are added smallest projected_perimeter_score first. The queue
// starts with both sides of every tunnel-weighted forest edge. Left of
// directed edge a→b, with a < b, is the side where (b-a) × (c-a) agrees
// with the edge normal. Accepting a triangle queues the unfilled side of
// each of its other edges that is not already under consideration. A
// rejected triangle is replaced by the next score on the same side. A side
// that already holds a triangle is closed. Every edge of a chosen triangle
// records the side that triangle occupies. A new candidate on an occupied
// side of any of its edges is skipped, and the next score is tried. A
// candidate that would be the third triangle on any of its edges is skipped
// the same way. So is one whose wedge at a vertex overlaps a triangle
// already incident there. Edges are projected onto the plane perpendicular
// to the vertex normal, and the wedges overlap when a non-shared edge of
// either triangle lies in the other's angle. The stored winding puts the
// edge normal on a newly chosen triangle, so the right triangle uses the
// opposite direction along the edge. After the queue is empty, each
// connected component is oriented from its first triangle. A neighbor whose
// cyclic order uses a shared edge in the same direction is reversed. That
// pass does not consult normals.
// Returns how many accepted triangles were proposed by an edge outside the
// forest.
struct Proposal {
    double score = 0.0;
    size_t corner = 0;
    size_t a = 0;
    size_t b = 0;
    int side = 0;
};

struct LaterProposal {
    // Smaller score first. Ties break by corner, then by the edge, then by side.
    bool operator()(const Proposal& left, const Proposal& right) const
    {
        if (left.score != right.score)
            return left.score > right.score;
        if (left.corner != right.corner)
            return left.corner > right.corner;
        if (left.a != right.a)
            return left.a > right.a;
        if (left.b != right.b)
            return left.b > right.b;
        return left.side > right.side;
    }
};

// One live proposal per open side. Slot 0 is the left side and slot 1 the right.
struct SideState {
    std::vector<std::pair<double, size_t>> ranked;
    size_t cursor = 0;
    bool closed = false;
    bool queued = false;
};

struct EdgeState {
    size_t a = 0;
    size_t b = 0;
    std::array<SideState, 2> side;
};

int side_slot(int side)
{
    return side > 0 ? 0 : 1;
}

struct TriangleFront {
    const std::vector<Vec3>& points;
    const std::vector<Vec3>& normals;
    const Geometry::AMGraph3D& graph;
    std::unordered_set<uint64_t> forest;

    std::unordered_map<TriKey, char, TriKeyHash> emitted;
    std::vector<std::array<size_t, 3>> faces;
    std::unordered_map<uint64_t, std::vector<size_t>> faces_on;
    std::unordered_map<uint64_t, unsigned char> side_taken;
    std::vector<std::vector<std::pair<size_t, size_t>>> vertex_tris;
    std::unordered_map<uint64_t, EdgeState> edges;
    std::priority_queue<Proposal, std::vector<Proposal>, LaterProposal> pending;
    int accepted_from_new_edges = 0;

    TriangleFront(const std::vector<Vec3>& points_in,
                  const std::vector<Vec3>& normals_in,
                  const Geometry::AMGraph3D& graph_in,
                  std::unordered_set<uint64_t> forest_in)
        : points(points_in)
        , normals(normals_in)
        , graph(graph_in)
        , forest(std::move(forest_in))
        , vertex_tris(points_in.size())
    {
    }

    // `a < b`. The higher index flips when the endpoint normals disagree.
    bool edge_frame(size_t a, size_t b, Vec3& edge_normal, double& normal_length) const
    {
        Vec3 normal_a = normals[a];
        Vec3 normal_b = normals[b];
        if (!finite_vector(normal_a) || !finite_vector(normal_b))
            return false;
        if (CGLA::dot(normal_a, normal_b) < 0.0)
            normal_b = -normal_b;
        edge_normal = normal_a + normal_b;
        normal_length = edge_normal.length();
        return normal_length > 1e-12 && std::isfinite(normal_length);
    }

    int edge_side(size_t u, size_t v, size_t corner) const
    {
        if (u > v)
            std::swap(u, v);
        if (corner == u || corner == v)
            return 0;
        Vec3 edge_normal;
        double normal_length = 0.0;
        if (!edge_frame(u, v, edge_normal, normal_length))
            return 0;
        const double area = triangle_area(points[u], points[v], points[corner]);
        if (!(area > 0.0) || !std::isfinite(area))
            return 0;
        const double side = CGLA::dot(
            CGLA::cross(points[v] - points[u], points[corner] - points[u]), edge_normal);
        if (!std::isfinite(side))
            return 0;
        const double cosine = side / (2.0 * area * normal_length);
        if (cosine > kSideEpsilon)
            return 1;
        if (cosine < -kSideEpsilon)
            return -1;
        return 0;
    }

    void mark_side(size_t u, size_t v, int side)
    {
        if (side == 0)
            return;
        side_taken[undirected_id(u, v)] |= side > 0 ? 1 : 2;
    }

    bool side_occupied(size_t u, size_t v, int side) const
    {
        if (side == 0)
            return false;
        const auto found = side_taken.find(undirected_id(u, v));
        if (found == side_taken.end())
            return false;
        const unsigned char bit = side > 0 ? 1 : 2;
        return (found->second & bit) != 0;
    }

    // The candidate occupies one side of edge (u, v). Reject it when that side
    // already holds a triangle.
    bool side_blocked(size_t u, size_t v, size_t other) const
    {
        return side_occupied(u, v, edge_side(u, v, other));
    }

    int incident_faces(size_t u, size_t v) const
    {
        const auto found = faces_on.find(undirected_id(u, v));
        if (found == faces_on.end())
            return 0;
        return static_cast<int>(found->second.size());
    }

    // Project the edges onto the plane perpendicular to the vertex normal.
    // The spoke is between the other two when that projection is a positive
    // combination of theirs. That is the wedge of at most 180 degrees, which
    // is the triangle's corner. Flipping the normal leaves the test unchanged.
    // An angle that collapses in the tangent plane does not block, and the
    // boundary of the angle is not inside.
    bool lies_between(size_t vertex, size_t end_a, size_t end_b, size_t spoke) const
    {
        const Vec3& normal = normals[vertex];
        if (!finite_vector(normal))
            return false;
        const double normal_length2 = CGLA::dot(normal, normal);
        if (!(normal_length2 > 1e-24) || !std::isfinite(normal_length2))
            return false;
        const Vec3 origin = points[vertex];
        auto tangent = [&](size_t end) {
            const Vec3 edge = points[end] - origin;
            return edge - normal * (CGLA::dot(edge, normal) / normal_length2);
        };
        const Vec3 edge_a = tangent(end_a);
        const Vec3 edge_b = tangent(end_b);
        const Vec3 edge_s = tangent(spoke);
        const double denom = CGLA::dot(CGLA::cross(edge_a, edge_b), normal);
        const double span = edge_a.length() * edge_b.length() * std::sqrt(normal_length2);
        if (!(span > 0.0) || !std::isfinite(span) || !std::isfinite(denom)
            || !(std::abs(denom) > kSideEpsilon * span))
            return false;
        const double alpha = CGLA::dot(CGLA::cross(edge_s, edge_b), normal) / denom;
        const double beta = CGLA::dot(CGLA::cross(edge_a, edge_s), normal) / denom;
        if (!std::isfinite(alpha) || !std::isfinite(beta))
            return false;
        return alpha > kSideEpsilon && beta > kSideEpsilon;
    }

    // The two sectors overlap when either triangle has a non-shared edge in
    // the other's wedge. A new triangle that contains an existing one has
    // both of its own edges outside that narrow wedge, so both directions
    // are required. A shared edge is the boundary and is allowed.
    bool vertex_blocked(size_t vertex, size_t spoke_u, size_t spoke_v) const
    {
        for (const auto& triangle : vertex_tris[vertex]) {
            const size_t end_a = triangle.first;
            const size_t end_b = triangle.second;
            if (spoke_u != end_a && spoke_u != end_b && lies_between(vertex, end_a, end_b, spoke_u))
                return true;
            if (spoke_v != end_a && spoke_v != end_b && lies_between(vertex, end_a, end_b, spoke_v))
                return true;
            if (end_a != spoke_u && end_a != spoke_v && lies_between(vertex, spoke_u, spoke_v, end_a))
                return true;
            if (end_b != spoke_u && end_b != spoke_v && lies_between(vertex, spoke_u, spoke_v, end_b))
                return true;
        }
        return false;
    }

    bool admissible(size_t a, size_t b, size_t corner) const
    {
        if (side_blocked(a, b, corner) || side_blocked(a, corner, b) || side_blocked(b, corner, a))
            return false;
        if (incident_faces(a, b) >= 2 || incident_faces(a, corner) >= 2
            || incident_faces(b, corner) >= 2)
            return false;
        return !vertex_blocked(a, b, corner) && !vertex_blocked(b, a, corner)
            && !vertex_blocked(corner, a, b);
    }

    bool emit(size_t i, size_t j, size_t k)
    {
        if (!emitted.emplace(tri_key(i, j, k), 1).second)
            return false;
        const size_t index = faces.size();
        faces.push_back({i, j, k});
        faces_on[undirected_id(i, j)].push_back(index);
        faces_on[undirected_id(j, k)].push_back(index);
        faces_on[undirected_id(k, i)].push_back(index);
        vertex_tris[i].emplace_back(j, k);
        vertex_tris[j].emplace_back(k, i);
        vertex_tris[k].emplace_back(i, j);
        mark_side(i, j, edge_side(i, j, k));
        mark_side(j, k, edge_side(j, k, i));
        mark_side(k, i, edge_side(k, i, j));
        return true;
    }

    // Corners that survive the tests that do not depend on the mesh, best first.
    void rank_side(EdgeState& edge, int side_sign, const Vec3& edge_normal, double normal_length)
    {
        SideState& side = edge.side[side_slot(side_sign)];
        const Vec3& point_a = points[edge.a];
        const Vec3& point_b = points[edge.b];
        const Vec3 along = point_b - point_a;
        for (Geometry::AMGraph::NodeID neighbor : graph.shared_neighbors(edge.a, edge.b)) {
            const size_t corner = neighbor;
            if (corner == edge.a || corner == edge.b)
                continue;
            const Vec3& point_c = points[corner];
            const double area = triangle_area(point_a, point_b, point_c);
            if (!(area > 0.0) || !std::isfinite(area))
                continue;
            const double signed_side = CGLA::dot(CGLA::cross(along, point_c - point_a), edge_normal);
            if (!std::isfinite(signed_side))
                continue;
            const double cosine = signed_side / (2.0 * area * normal_length);
            const bool on_side = side_sign > 0 ? cosine > kSideEpsilon : cosine < -kSideEpsilon;
            if (!on_side)
                continue;
            const double score = projected_perimeter_score(point_a, point_b, point_c, edge_normal);
            if (!std::isfinite(score))
                continue;
            side.ranked.emplace_back(score, corner);
        }
        std::sort(side.ranked.begin(), side.ranked.end());
    }

    // Push the best corner of this side that the mesh still allows. One entry
    // stays live until it is popped.
    void enqueue_next(uint64_t id, int side_sign)
    {
        const auto found = edges.find(id);
        if (found == edges.end())
            return;
        EdgeState& edge = found->second;
        SideState& side = edge.side[side_slot(side_sign)];
        if (side.closed || side.queued)
            return;
        if (side_occupied(edge.a, edge.b, side_sign)) {
            side.closed = true;
            return;
        }
        while (side.cursor < side.ranked.size()) {
            const size_t corner = side.ranked[side.cursor].second;
            if (emitted.find(tri_key(edge.a, edge.b, corner)) != emitted.end()) {
                mark_side(edge.a, edge.b, side_sign);
                side.closed = true;
                return;
            }
            if (admissible(edge.a, edge.b, corner)) {
                pending.push(Proposal{
                    side.ranked[side.cursor].first, corner, edge.a, edge.b, side_sign});
                side.queued = true;
                return;
            }
            ++side.cursor;
        }
        side.closed = true;
    }

    void seed()
    {
        for (const uint64_t id : forest) {
            EdgeState created;
            created.a = static_cast<size_t>(id >> 32);
            created.b = static_cast<size_t>(id & 0xffffffffu);
            Vec3 edge_normal;
            double normal_length = 0.0;
            if (!edge_frame(created.a, created.b, edge_normal, normal_length)) {
                created.side[0].closed = true;
                created.side[1].closed = true;
                edges.emplace(id, std::move(created));
                continue;
            }
            rank_side(created, 1, edge_normal, normal_length);
            rank_side(created, -1, edge_normal, normal_length);
            edges.emplace(id, std::move(created));
            enqueue_next(id, 1);
            enqueue_next(id, -1);
        }
    }

    // `other` is the third vertex of the triangle just added. Queue the side
    // that triangle does not already occupy. An edge already in the front is
    // left alone: its open side was queued when it entered.
    void open_created_edge(size_t u, size_t v, size_t other)
    {
        if (u == v)
            return;
        if (u > v)
            std::swap(u, v);
        const uint64_t id = undirected_id(u, v);
        if (edges.find(id) != edges.end())
            return;
        EdgeState created;
        created.a = u;
        created.b = v;
        Vec3 edge_normal;
        double normal_length = 0.0;
        const bool framed = edge_frame(u, v, edge_normal, normal_length);
        const int occupied = framed ? edge_side(u, v, other) : 0;
        if (!framed || occupied == 0) {
            created.side[0].closed = true;
            created.side[1].closed = true;
            edges.emplace(id, std::move(created));
            return;
        }
        created.side[side_slot(occupied)].closed = true;
        rank_side(created, -occupied, edge_normal, normal_length);
        edges.emplace(id, std::move(created));
        enqueue_next(id, -occupied);
    }

    // +1 when the cyclic order goes u then v, -1 when it goes v then u.
    static int edge_direction(const std::array<size_t, 3>& face, size_t u, size_t v)
    {
        for (int i = 0; i < 3; ++i) {
            const size_t from = face[static_cast<size_t>(i)];
            const size_t to = face[static_cast<size_t>((i + 1) % 3)];
            if (from == u && to == v)
                return 1;
            if (from == v && to == u)
                return -1;
        }
        return 0;
    }

    // The earliest triangle of a component keeps its order. The walk is the
    // recursive flood fill, stored on a stack so a long component does not
    // depend on the call stack.
    void orient()
    {
        const size_t count = faces.size();
        std::vector<char> seen(count, 0);
        std::vector<size_t> pending_faces;
        pending_faces.reserve(count);
        for (size_t seed = 0; seed < count; ++seed) {
            if (seen[seed])
                continue;
            seen[seed] = 1;
            pending_faces.push_back(seed);
            while (!pending_faces.empty()) {
                const size_t current = pending_faces.back();
                pending_faces.pop_back();
                const std::array<size_t, 3> face = faces[current];
                for (int edge = 0; edge < 3; ++edge) {
                    const size_t u = face[static_cast<size_t>(edge)];
                    const size_t v = face[static_cast<size_t>((edge + 1) % 3)];
                    const auto found = faces_on.find(undirected_id(u, v));
                    if (found == faces_on.end())
                        continue;
                    for (size_t neighbor : found->second) {
                        if (neighbor == current || seen[neighbor])
                            continue;
                        if (edge_direction(faces[neighbor], u, v) > 0)
                            std::swap(faces[neighbor][1], faces[neighbor][2]);
                        seen[neighbor] = 1;
                        pending_faces.push_back(neighbor);
                    }
                }
            }
        }
    }

    int grow()
    {
        seed();
        while (!pending.empty()) {
            const Proposal proposal = pending.top();
            pending.pop();
            const uint64_t id = undirected_id(proposal.a, proposal.b);
            const auto found = edges.find(id);
            if (found == edges.end())
                continue;
            SideState& side = found->second.side[side_slot(proposal.side)];
            side.queued = false;
            const bool have_triangle =
                emitted.find(tri_key(proposal.a, proposal.b, proposal.corner)) != emitted.end();
            if (side.closed || side_occupied(proposal.a, proposal.b, proposal.side) || have_triangle) {
                if (have_triangle || side_occupied(proposal.a, proposal.b, proposal.side))
                    mark_side(proposal.a, proposal.b, proposal.side);
                side.closed = true;
                continue;
            }
            if (!admissible(proposal.a, proposal.b, proposal.corner)) {
                ++side.cursor;
                enqueue_next(id, proposal.side);
                continue;
            }
            const bool added = proposal.side > 0
                ? emit(proposal.a, proposal.b, proposal.corner)
                : emit(proposal.a, proposal.corner, proposal.b);
            mark_side(proposal.a, proposal.b, proposal.side);
            side.closed = true;
            if (!added)
                continue;
            if (forest.count(id) == 0)
                ++accepted_from_new_edges;
            open_created_edge(proposal.a, proposal.corner, proposal.b);
            open_created_edge(proposal.b, proposal.corner, proposal.a);
        }
        return accepted_from_new_edges;
    }
};

int add_forest_triangles(const std::vector<Vec3>& points, const LRsROpts& opts, Manifold& result)
{
    const double spacing = median_neighbor_spacing(points);
    const double radius = (opts.radius_factor > 0.0) ? opts.radius_factor * spacing : 0.0;
    const Geometry::AMGraph3D graph = Geometry::delaunay_edges(points, radius);

    RSR::RSROpts normal_opts;
    normal_opts.dist = RSR::Distance::Euclidean;
    std::vector<Vec3> normals;
    RSR::point_cloud_normals(points, normal_opts, normals, RSR::NormalOrientation::None);
    if (normals.size() != points.size())
        normals.assign(points.size(), Vec3(0));

    const Geometry::GraphEdgeWeights weights = Geometry::tunnel_weight_graph(graph, normals);
    const Geometry::AMGraph3D forest = Geometry::minimum_spanning_forest(graph, weights.length);

    std::unordered_set<uint64_t> forest_ids;
    for (Geometry::AMGraph::NodeID node : forest.node_ids()) {
        for (const auto& entry : forest.edges(node)) {
            if (node < entry.first)
                forest_ids.insert(undirected_id(node, entry.first));
        }
    }

    TriangleFront front(points, normals, graph, std::move(forest_ids));
    const int accepted = front.grow();
    front.orient();
    for (const auto& face : front.faces) {
        const std::array<Vec3, 3> corners{
            points[face[0]], points[face[1]], points[face[2]]};
        result.add_face(corners);
    }
    stitch_mesh(result, 1e-8);
    result.cleanup();
    return accepted;
}

} // namespace

int reconstruct(const std::vector<Vec3>& points, const LRsROpts& opts, Manifold& result)
{
    result.clear();
    if (points.empty() || points.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        std::cout << "LRsR iterations: 0" << std::endl;
        return 0;
    }

    const int iterations = add_forest_triangles(points, opts, result);
    std::cout << "LRsR iterations: " << iterations << std::endl;
    return iterations;
}

} // namespace HMesh::LRsR
