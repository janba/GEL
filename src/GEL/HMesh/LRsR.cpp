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
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace HMesh::LRsR {
namespace {

using NodeID = Geometry::AMGraph::NodeID;
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

/// One Delaunay link, kept whole. A closed link is an umbrella. An open link
/// is the partial disk those triangles cover.
struct IncidentDisk {
    bool closed = false;
    double area = 0;
    std::vector<int> neighbors;
    std::vector<std::pair<size_t, size_t>> ends;
    std::vector<char> forced;
};

bool disk_has_node(const IncidentDisk& disk, NodeID node)
{
    for (int neighbor : disk.neighbors)
        if (static_cast<NodeID>(neighbor) == node)
            return true;
    return false;
}

bool disk_has_triangle(const IncidentDisk& disk, size_t vertex, const TriKey& key)
{
    for (const auto& end : disk.ends)
        if (tri_key(vertex, end.first, end.second) == key)
            return true;
    return false;
}

/// Full umbrellas come first. Area then chooses, and a lower index breaks a tie.
bool disk_precedes(const std::vector<IncidentDisk>& disks, int i, int j)
{
    if (disks[static_cast<size_t>(i)].closed != disks[static_cast<size_t>(j)].closed)
        return disks[static_cast<size_t>(i)].closed;
    const double ai = disks[static_cast<size_t>(i)].area;
    const double aj = disks[static_cast<size_t>(j)].area;
    if (ai < aj)
        return true;
    if (aj < ai)
        return false;
    return i < j;
}

/// Keep a whole link. A closed link that already contains every incident tree
/// neighbor and every earlier triangle is the disk, and the smallest such
/// area wins when several links do. Otherwise the vertex still receives one
/// minimum-area umbrella when any link is closed, and further whole links
/// until the remaining tree edges and earlier triangles are covered. With no
/// closed link, the same rule keeps a partial disk.
std::vector<char> disks_to_keep(const std::vector<IncidentDisk>& disks,
                                const std::vector<NodeID>& tree_neighbors,
                                size_t vertex)
{
    const int n = static_cast<int>(disks.size());
    std::vector<char> keep(static_cast<size_t>(n), 0);
    if (n == 0)
        return keep;

    std::vector<NodeID> pending_nodes;
    for (NodeID node : tree_neighbors) {
        for (const IncidentDisk& disk : disks) {
            if (disk_has_node(disk, node)) {
                pending_nodes.push_back(node);
                break;
            }
        }
    }
    std::vector<TriKey> pending_tris;
    for (const IncidentDisk& disk : disks) {
        for (size_t t = 0; t < disk.ends.size(); ++t) {
            if (!disk.forced[t])
                continue;
            const TriKey key = tri_key(vertex, disk.ends[t].first, disk.ends[t].second);
            if (std::find(pending_tris.begin(), pending_tris.end(), key) == pending_tris.end())
                pending_tris.push_back(key);
        }
    }

    const auto consume = [&](const IncidentDisk& disk) {
        std::erase_if(pending_nodes, [&](NodeID node) { return disk_has_node(disk, node); });
        std::erase_if(pending_tris, [&](const TriKey& key) {
            return disk_has_triangle(disk, vertex, key);
        });
    };
    const auto covers = [&](const IncidentDisk& disk) {
        for (NodeID node : pending_nodes)
            if (disk_has_node(disk, node))
                return true;
        for (const TriKey& key : pending_tris)
            if (disk_has_triangle(disk, vertex, key))
                return true;
        return false;
    };
    const auto covers_every_constraint = [&](const IncidentDisk& disk) {
        for (NodeID node : pending_nodes)
            if (!disk_has_node(disk, node))
                return false;
        for (const TriKey& key : pending_tris)
            if (!disk_has_triangle(disk, vertex, key))
                return false;
        return true;
    };

    int best_full = -1;
    for (int i = 0; i < n; ++i) {
        if (!disks[static_cast<size_t>(i)].closed || !covers_every_constraint(disks[static_cast<size_t>(i)]))
            continue;
        if (best_full < 0 || disk_precedes(disks, i, best_full))
            best_full = i;
    }
    if (best_full >= 0) {
        keep[static_cast<size_t>(best_full)] = 1;
        return keep;
    }

    int any_full = -1;
    for (int i = 0; i < n; ++i) {
        if (disks[static_cast<size_t>(i)].closed && (any_full < 0 || disk_precedes(disks, i, any_full)))
            any_full = i;
    }
    if (any_full >= 0) {
        keep[static_cast<size_t>(any_full)] = 1;
        consume(disks[static_cast<size_t>(any_full)]);
    }
    for (int step = 0; step < n; ++step) {
        int best = -1;
        for (int i = 0; i < n; ++i) {
            if (keep[static_cast<size_t>(i)] || !covers(disks[static_cast<size_t>(i)]))
                continue;
            if (best < 0 || disk_precedes(disks, i, best))
                best = i;
        }
        if (best < 0)
            break;
        keep[static_cast<size_t>(best)] = 1;
        consume(disks[static_cast<size_t>(best)]);
    }
    if (std::find(keep.begin(), keep.end(), 1) == keep.end()) {
        int best = 0;
        for (int i = 1; i < n; ++i)
            if (disk_precedes(disks, i, best))
                best = i;
        keep[static_cast<size_t>(best)] = 1;
    }
    return keep;
}

bool next_is(const std::array<size_t, 3>& face, size_t vertex, size_t next)
{
    for (int i = 0; i < 3; ++i)
        if (face[static_cast<size_t>(i)] == vertex)
            return face[static_cast<size_t>((i + 1) % 3)] == next;
    return false;
}

struct PairHash {
    size_t operator()(const std::pair<size_t, size_t>& edge) const noexcept
    {
        size_t h = edge.first;
        h ^= edge.second + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

/// The proposed umbrellas can stack three triangles on one edge, or wind an
/// edge twice in the same direction. Those are extra triangles. Drop the
/// larger one until every edge has one winding and at most two triangles.
void prune_extra_triangles(std::vector<std::array<size_t, 3>>& faces,
                           const std::vector<Vec3>& points)
{
    const int n = static_cast<int>(faces.size());
    std::vector<char> alive(static_cast<size_t>(n), 1);
    std::vector<double> area(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        const auto& face = faces[static_cast<size_t>(i)];
        area[static_cast<size_t>(i)] = triangle_area(points[face[0]], points[face[1]], points[face[2]]);
    }

    // True when face i is the one to drop: larger area, then later index.
    const auto worse = [&](int i, int j) {
        const double ai = area[static_cast<size_t>(i)];
        const double aj = area[static_cast<size_t>(j)];
        if (aj < ai)
            return true;
        if (ai < aj)
            return false;
        return i > j;
    };

    for (int step = 0; step < n; ++step) {
        std::unordered_map<std::pair<size_t, size_t>, std::vector<int>, PairHash> undirected;
        std::unordered_map<std::pair<size_t, size_t>, std::vector<int>, PairHash> directed;
        for (int i = 0; i < n; ++i) {
            if (!alive[static_cast<size_t>(i)])
                continue;
            const auto& face = faces[static_cast<size_t>(i)];
            for (int e = 0; e < 3; ++e) {
                size_t u = face[static_cast<size_t>(e)];
                size_t v = face[static_cast<size_t>((e + 1) % 3)];
                directed[{u, v}].push_back(i);
                if (u > v)
                    std::swap(u, v);
                undirected[{u, v}].push_back(i);
            }
        }
        int victim = -1;
        for (const auto& entry : undirected) {
            if (entry.second.size() < 3)
                continue;
            for (int idx : entry.second)
                if (victim < 0 || worse(idx, victim))
                    victim = idx;
        }
        for (const auto& entry : directed) {
            if (entry.second.size() < 2)
                continue;
            for (int idx : entry.second)
                if (victim < 0 || worse(idx, victim))
                    victim = idx;
        }
        if (victim < 0)
            break;
        alive[static_cast<size_t>(victim)] = 0;
    }

    std::vector<std::array<size_t, 3>> kept;
    kept.reserve(faces.size());
    for (int i = 0; i < n; ++i)
        if (alive[static_cast<size_t>(i)])
            kept.push_back(faces[static_cast<size_t>(i)]);
    faces.swap(kept);
}

} // namespace

void reconstruct(const std::vector<Vec3>& points, const LRsROpts& opts, Manifold& result)
{
    result.clear();
    const size_t n = points.size();
    if (n == 0)
        return;

    const double spacing = median_neighbor_spacing(points);
    const double radius = (opts.radius_factor > 0.0) ? opts.radius_factor * spacing : 0.0;
    const std::vector<Geometry::DelaunayStar> stars = Geometry::delaunay_stars(points, radius);

    Geometry::AMGraph3D delaunay;
    for (size_t i = 0; i < n; ++i)
        delaunay.add_node(points[i]);
    if (stars.size() == n) {
        for (NodeID node = 0; node < n; ++node) {
            for (int neighbor : stars[node].neighbors)
                delaunay.connect_nodes(node, static_cast<NodeID>(neighbor));
        }
    }

    std::vector<Vec3> normals;
    RSR::RSROpts normal_opts;
    normal_opts.dist = RSR::Distance::Euclidean;
    RSR::point_cloud_normals(points, normal_opts, normals, RSR::NormalOrientation::None);
    if (normals.size() != n)
        normals.assign(n, Vec3(0, 0, 1));

    const Geometry::GraphEdgeWeights weights = Geometry::tunnel_weight_graph(delaunay, normals);
    const Geometry::AMGraph3D forest = Geometry::minimum_spanning_forest(delaunay, weights.length);

    std::vector<std::vector<NodeID>> tree_nbr(n);
    if (forest.no_nodes() == n) {
        for (NodeID a = 0; a < n; ++a) {
            for (NodeID b : forest.neighbors(a)) {
                if (a < b)
                    tree_nbr[a].push_back(b), tree_nbr[b].push_back(a);
            }
        }
    }
    for (auto& adj : tree_nbr)
        std::sort(adj.begin(), adj.end());

    std::unordered_map<TriKey, std::array<size_t, 3>, TriKeyHash> chosen;
    std::vector<std::array<size_t, 3>> faces;
    std::vector<char> seen(n, 0);
    std::queue<NodeID> queue;

    for (NodeID root = 0; root < n; ++root) {
        if (seen[root])
            continue;
        seen[root] = 1;
        queue.push(root);
        while (!queue.empty()) {
            const NodeID vertex = queue.front();
            queue.pop();

            if (vertex < stars.size()) {
                std::vector<IncidentDisk> disks;
                for (const Geometry::DelaunayLink& comp : stars[vertex].links) {
                    const int k = static_cast<int>(comp.neighbors.size());
                    const bool closed = comp.closed && k >= 3;
                    const int m = closed ? k : k - 1;
                    if (m < 1)
                        continue;

                    IncidentDisk disk;
                    disk.closed = closed;
                    disk.neighbors = comp.neighbors;
                    disk.ends.resize(static_cast<size_t>(m));
                    disk.forced.resize(static_cast<size_t>(m), 0);
                    for (int t = 0; t < m; ++t) {
                        const int i1 = closed ? (t + 1) % k : t + 1;
                        const size_t a = static_cast<size_t>(comp.neighbors[static_cast<size_t>(t)]);
                        const size_t b = static_cast<size_t>(comp.neighbors[static_cast<size_t>(i1)]);
                        disk.ends[static_cast<size_t>(t)] = {a, b};
                        disk.area += triangle_area(points[vertex], points[a], points[b]);
                        if (chosen.contains(tri_key(vertex, a, b)))
                            disk.forced[static_cast<size_t>(t)] = 1;
                    }
                    disks.push_back(std::move(disk));
                }

                const std::vector<char> keep = disks_to_keep(disks, tree_nbr[vertex], vertex);
                for (size_t di = 0; di < disks.size(); ++di) {
                    if (!keep[di])
                        continue;
                    const IncidentDisk& disk = disks[di];
                    bool flip = false;
                    for (size_t t = 0; t < disk.ends.size(); ++t) {
                        if (!disk.forced[t])
                            continue;
                        const auto [a, b] = disk.ends[t];
                        const auto it = chosen.find(tri_key(vertex, a, b));
                        // The stored successor of this vertex has to be `a`, the
                        // first corner in the link order. Otherwise the new fan
                        // is wound the other way so it matches that face.
                        if (it != chosen.end() && !next_is(it->second, vertex, a))
                            flip = true;
                        break;
                    }
                    for (const auto& end : disk.ends) {
                        const auto [a, b] = end;
                        const std::array<size_t, 3> face = flip
                            ? std::array<size_t, 3>{vertex, b, a}
                            : std::array<size_t, 3>{vertex, a, b};
                        const TriKey key = tri_key(face[0], face[1], face[2]);
                        if (chosen.emplace(key, face).second)
                            faces.push_back(face);
                    }
                }
            }

            for (NodeID neighbor : tree_nbr[vertex]) {
                if (neighbor < n && !seen[neighbor]) {
                    seen[neighbor] = 1;
                    queue.push(neighbor);
                }
            }
        }
    }

    prune_extra_triangles(faces, points);
    for (const std::array<size_t, 3>& face : faces) {
        const std::array<Vec3, 3> corners{
            points[face[0]], points[face[1]], points[face[2]]};
        result.add_face(corners);
    }
    // Copies of one input point land on the same coordinates. A positive
    // radius clusters those copies; the strict ball test excludes distance 0
    // when the radius is 0.
    stitch_mesh(result, 1e-8);
    result.cleanup();
}

} // namespace HMesh::LRsR
