/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/GGSR.h>

#include <GEL/Geometry/Delaunay3D.h>
#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/cleanup.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace HMesh::GGSR {
namespace {

using Vec3 = CGLA::Vec3d;
using NodeID = Geometry::AMGraph::NodeID;

constexpr double kDirEps = 1e-24;
/// A turn this close to straight does not count as progress toward closing.
constexpr double kTurnEps = 1e-9;
/// One corner belongs to one cycle when this is set. Temporarily off so a
/// later walk can still cross a corner that an earlier face already used.
constexpr bool kExclusiveCorners = false;

struct EdgeKey {
    NodeID lo = 0;
    NodeID hi = 0;
    bool operator==(const EdgeKey& other) const
    {
        return lo == other.lo && hi == other.hi;
    }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& edge) const noexcept
    {
        size_t h = edge.lo;
        h ^= edge.hi + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

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

bool walks(const std::vector<NodeID>& cycle, NodeID from, NodeID to)
{
    const size_t n = cycle.size();
    for (size_t i = 0; i < n; ++i)
        if (cycle[i] == from && cycle[(i + 1) % n] == to)
            return true;
    return false;
}

EdgeKey edge_key(NodeID a, NodeID b)
{
    if (a > b)
        std::swap(a, b);
    return {a, b};
}

/// The two neighbors that form a corner. Stored with `a < b`.
struct Corner {
    NodeID a = 0;
    NodeID b = 0;
};

/// Every unordered pair of incident edges with a defined direction. A pair
/// is stored with `a < b`.
std::vector<std::vector<Corner>> corners_at(const Geometry::AMGraph3D& graph)
{
    const size_t n = graph.no_nodes();
    std::vector<std::vector<Corner>> corners(n);
    std::vector<NodeID> ids;
    for (NodeID vertex = 0; vertex < n; ++vertex) {
        ids.clear();
        for (const NodeID neighbor : graph.neighbors_lazy(vertex)) {
            const Vec3 dir = graph.pos[neighbor] - graph.pos[vertex];
            const double len2 = CGLA::sqr_length(dir);
            if (!(len2 > kDirEps) || !std::isfinite(len2))
                continue;
            ids.push_back(neighbor);
        }
        if (ids.size() < 2)
            continue;
        auto& local = corners[vertex];
        local.reserve(ids.size() * (ids.size() - 1) / 2);
        for (size_t i = 0; i < ids.size(); ++i) {
            for (size_t j = i + 1; j < ids.size(); ++j)
                local.push_back(Corner{ids[i], ids[j]});
        }
    }
    return corners;
}

bool on_path(const std::vector<NodeID>& path, NodeID id)
{
    for (const NodeID have : path)
        if (have == id)
            return true;
    return false;
}

bool corner_taken(const std::vector<std::vector<Corner>>& used, NodeID vertex, NodeID a, NodeID b)
{
    if (!kExclusiveCorners)
        return false;
    if (a > b)
        std::swap(a, b);
    for (const Corner& corner : used[vertex])
        if (corner.a == a && corner.b == b)
            return true;
    return false;
}

void claim_corner(std::vector<std::vector<Corner>>& used, NodeID vertex, NodeID a, NodeID b)
{
    if (!kExclusiveCorners)
        return;
    if (corner_taken(used, vertex, a, b))
        return;
    if (a > b)
        std::swap(a, b);
    used[vertex].push_back(Corner{a, b});
}

void claim_cycle(std::vector<std::vector<Corner>>& used, const std::vector<NodeID>& cycle)
{
    const size_t n = cycle.size();
    for (size_t i = 0; i < n; ++i)
        claim_corner(used, cycle[i], cycle[(i + n - 1) % n], cycle[(i + 1) % n]);
}

/// Walk from `start` along `depart`. The other corner edge fixes the normal
/// so the corner lies to the left of travel. Each later step looks at every
/// other incident edge. A positive turn back to an earlier vertex of this
/// walk closes that shorter cycle.
/// Otherwise a positive turn whose far vertex neighbors `start` is taken,
/// because the cycle then closes on the next step. The shortest such cycle
/// wins. Otherwise the greatest signed turn wins, and an equal turn keeps
/// the smaller neighbor id. When `kExclusiveCorners` is set, a chosen corner
/// that is already claimed ends the walk. The accepted corner then supplies
/// the normal, flipped into the same hemisphere.
std::vector<NodeID> trace_corner(const Geometry::AMGraph3D& graph,
                                 const std::vector<std::vector<Corner>>& used,
                                 NodeID start,
                                 NodeID depart,
                                 NodeID other,
                                 int max_edges)
{
    const Vec3 from_depart = graph.pos[depart] - graph.pos[start];
    const Vec3 from_other = graph.pos[other] - graph.pos[start];
    Vec3 normal = CGLA::cross(from_depart, from_other);
    const double normal_len2 = CGLA::sqr_length(normal);
    if (!(normal_len2 > kDirEps) || !std::isfinite(normal_len2))
        return {};
    normal /= std::sqrt(normal_len2);

    std::vector<NodeID> path;
    path.reserve(static_cast<size_t>(max_edges));
    path.push_back(start);
    NodeID prev = start;
    NodeID cur = depart;
    int edges = 0;
    const int no_close = max_edges + 1;
    while (edges < max_edges) {
        ++edges;
        if (cur == start) {
            if (edges < 3 || path.size() < 2)
                return {};
            if (corner_taken(used, start, path[1], prev))
                return {};
            return std::move(path);
        }
        if (on_path(path, cur))
            return {};
        path.push_back(cur);

        const Vec3 inbound = graph.pos[cur] - graph.pos[prev];
        const Vec3 in_plane = inbound - normal * CGLA::dot(inbound, normal);
        const double in_len2 = CGLA::sqr_length(in_plane);
        if (!(in_len2 > kDirEps) || !std::isfinite(in_len2))
            return {};

        bool found = false;
        int best_close = no_close;
        double best_score = 0.0;
        NodeID best_next = 0;
        Vec3 best_cross(0);
        for (const NodeID next : graph.neighbors_lazy(cur)) {
            if (next == prev)
                continue;
            const Vec3 outbound = graph.pos[next] - graph.pos[cur];
            const Vec3 out_plane = outbound - normal * CGLA::dot(outbound, normal);
            const double out_len2 = CGLA::sqr_length(out_plane);
            if (!(out_len2 > kDirEps) || !std::isfinite(out_len2))
                continue;
            const double score = std::atan2(CGLA::dot(normal, CGLA::cross(in_plane, out_plane)),
                                            CGLA::dot(in_plane, out_plane));
            if (!std::isfinite(score))
                continue;
            int close = no_close;
            const auto earlier = std::find(path.begin(), path.end(), next);
            if (earlier != path.end()) {
                const int back = static_cast<int>(path.end() - earlier);
                // Repeating a vertex cannot extend the walk. It is useful only
                // when the turn closes the part of the walk from that vertex.
                if (!(score > kTurnEps) || back < 3 || back > max_edges ||
                    corner_taken(used, *earlier, *(earlier + 1), cur))
                    continue;
                close = back;
            } else if (score > kTurnEps) {
                const int close_next = edges + 2;
                if (close_next >= 3 && close_next <= max_edges && graph.edges(next).contains(start))
                    close = close_next;
            }
            if (!found || close < best_close ||
                (close == best_close &&
                 (score > best_score || (score == best_score && next < best_next)))) {
                found = true;
                best_close = close;
                best_score = score;
                best_next = next;
                best_cross = CGLA::cross(graph.pos[prev] - graph.pos[cur], outbound);
            }
        }
        if (!found || corner_taken(used, cur, prev, best_next))
            return {};
        const auto closing = std::find(path.begin(), path.end(), best_next);
        if (closing != path.end())
            return {closing, path.end()};

        const double cross_len2 = CGLA::sqr_length(best_cross);
        if (cross_len2 > kDirEps && std::isfinite(cross_len2)) {
            if (CGLA::dot(best_cross, normal) < 0.0)
                best_cross = -best_cross;
            normal = best_cross / std::sqrt(cross_len2);
        }
        prev = cur;
        cur = best_next;
    }
    return {};
}

void canonicalize(std::vector<NodeID>& cycle)
{
    const auto least = std::min_element(cycle.begin(), cycle.end());
    std::rotate(cycle.begin(), least, cycle.end());
    if (cycle.size() >= 2 && cycle[1] > cycle.back())
        std::reverse(cycle.begin() + 1, cycle.end());
}

std::vector<std::vector<NodeID>> trace_faces(const Geometry::AMGraph3D& graph, int max_edges)
{
    std::vector<std::vector<NodeID>> faces;
    if (max_edges < 3 || graph.no_nodes() == 0)
        return faces;
    int limit = max_edges;
    if (graph.no_nodes() < static_cast<size_t>(std::numeric_limits<int>::max()))
        limit = std::min(limit, static_cast<int>(graph.no_nodes()));
    const std::vector<std::vector<Corner>> corners = corners_at(graph);
    const size_t n = graph.no_nodes();
    std::vector<std::vector<Corner>> used(n);
    for (NodeID vertex = 0; vertex < n; ++vertex) {
        for (const Corner& corner : corners[vertex]) {
            if (corner_taken(used, vertex, corner.a, corner.b))
                continue;
            for (int side = 0; side < 2; ++side) {
                if (corner_taken(used, vertex, corner.a, corner.b))
                    break;
                const NodeID depart = side == 0 ? corner.a : corner.b;
                const NodeID other = side == 0 ? corner.b : corner.a;
                std::vector<NodeID> cycle = trace_corner(graph, used, vertex, depart, other, limit);
                if (cycle.size() < 3)
                    continue;
                claim_cycle(used, cycle);
                canonicalize(cycle);
                faces.push_back(std::move(cycle));
            }
        }
    }
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    return faces;
}

void orient_shared_edges(std::vector<std::vector<NodeID>>& faces)
{
    std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash> incident;
    for (size_t face = 0; face < faces.size(); ++face) {
        const std::vector<NodeID>& cycle = faces[face];
        const size_t n = cycle.size();
        for (size_t i = 0; i < n; ++i)
            incident[edge_key(cycle[i], cycle[(i + 1) % n])].push_back(face);
    }

    std::vector<char> seen(faces.size(), 0);
    std::queue<size_t> pending;
    for (size_t start = 0; start < faces.size(); ++start) {
        if (seen[start])
            continue;
        seen[start] = 1;
        pending.push(start);
        while (!pending.empty()) {
            const size_t face = pending.front();
            pending.pop();
            const std::vector<NodeID>& cycle = faces[face];
            const size_t n = cycle.size();
            for (size_t i = 0; i < n; ++i) {
                const NodeID from = cycle[i];
                const NodeID to = cycle[(i + 1) % n];
                const auto found = incident.find(edge_key(from, to));
                if (found == incident.end())
                    continue;
                for (const size_t other : found->second) {
                    if (seen[other])
                        continue;
                    if (walks(faces[other], from, to))
                        std::reverse(faces[other].begin(), faces[other].end());
                    seen[other] = 1;
                    pending.push(other);
                }
            }
        }
    }
}

} // namespace

int reconstruct(const std::vector<Vec3>& points, const GGSROpts& opts, Manifold& result)
{
    result.clear();
    if (points.empty() || points.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        std::cout << "GGSR mst edges: 0" << std::endl;
        std::cout << "GGSR faces: 0" << std::endl;
        return 0;
    }

    const double spacing = median_neighbor_spacing(points);
    const double radius = (opts.radius_factor > 0.0) ? opts.radius_factor * spacing : 0.0;
    const Geometry::AMGraph3D graph = Geometry::gabriel_edges(points, radius);

    RSR::RSROpts normal_opts;
    normal_opts.dist = RSR::Distance::Euclidean;
    std::vector<Vec3> normals;
    RSR::point_cloud_normals(points, normal_opts, normals, RSR::NormalOrientation::None);
    if (normals.size() != points.size())
        normals.assign(points.size(), Vec3(0));

    const Geometry::GraphEdgeWeights weights = Geometry::tunnel_weight_graph(graph, normals);
    const Geometry::AMGraph3D forest = Geometry::minimum_spanning_forest(graph, weights.length);

    std::vector<std::vector<NodeID>> faces = trace_faces(graph, opts.max_edges);
    orient_shared_edges(faces);

    for (const std::vector<NodeID>& face : faces) {
        std::vector<Vec3> corners;
        corners.reserve(face.size());
        for (const NodeID id : face)
            corners.push_back(points[id]);
        result.add_face(corners);
    }
    stitch_mesh(result, 1e-8);
    result.cleanup();

    const int count = static_cast<int>(faces.size());
    std::cout << "GGSR mst edges: " << forest.no_edges() << std::endl;
    std::cout << "GGSR faces: " << count << std::endl;
    return count;
}

} // namespace HMesh::GGSR
