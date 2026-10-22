/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/GGSR.h>

#include <GEL/Geometry/Delaunay3D.h>
#include <GEL/Geometry/LocalDijkstraFaces.h>
#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/cleanup.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

struct FaceMeta {
    NodeID start = 0;
    Vec3 bary;
    double perimeter = 0;
    size_t edge_count = 0;
};

struct Projected {
    Vec3 vector;
    bool ok = false;
};

bool strictly_larger(const FaceMeta& lhs, const FaceMeta& rhs)
{
    if (lhs.edge_count != rhs.edge_count)
        return lhs.edge_count > rhs.edge_count;
    return lhs.perimeter > rhs.perimeter;
}

Projected project_from_start(NodeID start, const Vec3& bary, const std::vector<Vec3>& points, const std::vector<Vec3>& normals)
{
    Projected projected;
    if (start >= points.size() || start >= normals.size())
        return projected;
    const Vec3& normal = normals[start];
    const double scale = CGLA::sqr_length(normal);
    if (!(scale > 0.0))
        return projected;
    const Vec3 offset = bary - points[start];
    projected.vector = offset - normal * (CGLA::dot(offset, normal) / scale);
    projected.ok = true;
    return projected;
}

/// Drop the strictly larger face when two faces share more than three edges.
/// The same happens when they share exactly two and the projected vectors from
/// their starting vertices to their barycenters have a positive dot product.
/// A missing normal leaves that projection undefined, so the pair is kept.
/// Shared edges are counted in arrays indexed by face. Pairs are not stored.
std::vector<std::vector<NodeID>> discard_larger_faces(std::vector<Geometry::DijkstraFace>& found,
                                                      const std::vector<Vec3>& points,
                                                      const std::vector<Vec3>& normals)
{
    const size_t face_count = found.size();
    std::vector<FaceMeta> meta(face_count);
    std::unordered_map<EdgeKey, uint32_t, EdgeKeyHash> edge_of;
    std::vector<std::vector<uint32_t>> users;
    std::vector<std::vector<uint32_t>> face_edges(face_count);

    for (size_t index = 0; index < face_count; ++index) {
        const std::vector<NodeID>& vertices = found[index].vertices;
        const size_t count = vertices.size();
        FaceMeta& item = meta[index];
        item.start = found[index].start;
        item.edge_count = count;
        Vec3 sum(0);
        auto& edges = face_edges[index];
        edges.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            const NodeID from = vertices[i];
            const NodeID to = vertices[(i + 1) % count];
            sum += points[from];
            item.perimeter += CGLA::length(points[to] - points[from]);
            const EdgeKey key = edge_key(from, to);
            const auto inserted = edge_of.emplace(key, static_cast<uint32_t>(edge_of.size()));
            const uint32_t edge = inserted.first->second;
            if (inserted.second)
                users.emplace_back();
            bool seen = false;
            for (const uint32_t have : edges) {
                if (have == edge) {
                    seen = true;
                    break;
                }
            }
            if (seen)
                continue;
            edges.push_back(edge);
            users[edge].push_back(static_cast<uint32_t>(index));
        }
        item.bary = count > 0 ? sum / static_cast<double>(count) : Vec3(0);
    }

    std::vector<char> drop(face_count, 0);
    if (face_count <= std::numeric_limits<uint32_t>::max()) {
        std::vector<unsigned> stamp(face_count, 0);
        std::vector<int> share(face_count, 0);
        std::vector<uint32_t> touched;
        touched.reserve(256);
        unsigned epoch = 0;
        for (uint32_t index = 0; index < static_cast<uint32_t>(face_count); ++index) {
            if (++epoch == 0) {
                std::fill(stamp.begin(), stamp.end(), 0);
                epoch = 1;
            }
            touched.clear();
            for (const uint32_t edge : face_edges[index]) {
                for (const uint32_t other : users[edge]) {
                    if (other == index)
                        continue;
                    if (stamp[other] != epoch) {
                        stamp[other] = epoch;
                        share[other] = 1;
                        touched.push_back(other);
                    } else
                        ++share[other];
                }
            }
            const FaceMeta& lhs = meta[index];
            for (const uint32_t other : touched) {
                if (other < index)
                    continue;
                const int count = share[other];
                bool reject = count > 3;
                if (count == 2) {
                    const Projected left = project_from_start(lhs.start, lhs.bary, points, normals);
                    const Projected right = project_from_start(meta[other].start, meta[other].bary, points, normals);
                    reject = left.ok && right.ok && CGLA::dot(left.vector, right.vector) > 0.0;
                }
                if (!reject)
                    continue;
                if (strictly_larger(lhs, meta[other]))
                    drop[index] = 1;
                else if (strictly_larger(meta[other], lhs))
                    drop[other] = 1;
            }
        }
    }

    std::vector<std::vector<NodeID>> kept;
    kept.reserve(face_count);
    for (size_t index = 0; index < face_count; ++index) {
        if (drop[index])
            continue;
        kept.push_back(std::move(found[index].vertices));
    }
    return kept;
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

    std::vector<Geometry::DijkstraFace> found = Geometry::local_dijkstra_faces(graph, opts.max_edges);
    std::vector<std::vector<NodeID>> faces = discard_larger_faces(found, points, normals);
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
