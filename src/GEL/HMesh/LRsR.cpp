/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/LRsR.h>

#include <GEL/Geometry/Delaunay3D.h>
#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/normal.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/cleanup.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
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
    bool operator<(const TriKey& other) const
    {
        if (a != other.a) return a < other.a;
        if (b != other.b) return b < other.b;
        return c < other.c;
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
    return perimeter * area / projected_area;
}

// Six nearest neighbors. One of them may lie on another sheet, so the normal
// is the most planar leave-one-out fit. A bridge triangle leans more than
// 60 degrees off that normal.
constexpr int kNormalNeighbors = 6;
constexpr double kNormalAgreement = 0.5;

struct PlanarFit {
    Vec3 normal = Vec3(0);
    double flatness = 1e9;
};

PlanarFit fit_points(const std::vector<Vec3>& samples)
{
    PlanarFit fit;
    if (samples.size() < 3)
        return fit;
    Vec3 centroid(0);
    for (const Vec3& sample : samples)
        centroid += sample;
    centroid /= static_cast<double>(samples.size());
    CGLA::Mat3x3d covariance(0);
    for (const Vec3& sample : samples) {
        const Vec3 diff = sample - centroid;
        covariance += CGLA::outer_product(diff, diff);
    }
    CGLA::Mat3x3d eigenvectors(0);
    std::array<double, 3> eigenvalues{};
    Geometry::eigen3x3(covariance, eigenvectors, eigenvalues);
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int a, int b) { return eigenvalues[a] < eigenvalues[b]; });
    fit.normal = eigenvectors[order[0]];
    const double middle = std::max(eigenvalues[order[1]], 1e-18);
    fit.flatness = std::max(0.0, eigenvalues[order[0]]) / middle;
    return fit;
}

std::vector<Vec3> local_normals(const std::vector<Vec3>& points)
{
    std::vector<Vec3> normals(points.size(), Vec3(0));
    if (points.size() < 4)
        return normals;
    Geometry::Tree tree;
    for (size_t i = 0; i < points.size(); ++i) {
        const Vec3& point = points[i];
        if (std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]))
            tree.insert(point, i);
    }
    tree.build();

    std::vector<Geometry::Tree::IdRecord> hits;
    std::vector<Vec3> neighbors;
    neighbors.reserve(static_cast<size_t>(kNormalNeighbors));
    for (size_t i = 0; i < points.size(); ++i) {
        hits.clear();
        tree.m_closest_ids(static_cast<unsigned>(kNormalNeighbors + 1), points[i],
                           std::numeric_limits<double>::infinity(), hits);
        std::sort(hits.begin(), hits.end(),
                  [](const Geometry::Tree::IdRecord& a, const Geometry::Tree::IdRecord& b) {
                      return a.d < b.d;
                  });
        neighbors.clear();
        for (const Geometry::Tree::IdRecord& hit : hits) {
            if (hit.v == i || !(hit.d > 1e-24) || !std::isfinite(hit.d))
                continue;
            if (static_cast<int>(neighbors.size()) >= kNormalNeighbors)
                break;
            neighbors.push_back(points[static_cast<size_t>(hit.v)]);
        }
        if (neighbors.size() < 3)
            continue;
        PlanarFit best = fit_points(neighbors);
        if (neighbors.size() >= 5) {
            std::vector<Vec3> subset;
            subset.reserve(neighbors.size() - 1);
            for (size_t drop = 0; drop < neighbors.size(); ++drop) {
                subset.clear();
                for (size_t j = 0; j < neighbors.size(); ++j)
                    if (j != drop)
                        subset.push_back(neighbors[j]);
                const PlanarFit fit = fit_points(subset);
                if (fit.flatness < best.flatness)
                    best = fit;
            }
        }
        const Vec3& normal = best.normal;
        if (!std::isfinite(normal[0]) || !std::isfinite(normal[1]) || !std::isfinite(normal[2]))
            continue;
        if (normal.length() < 0.5)
            continue;
        normals[static_cast<size_t>(i)] = normal;
    }
    return normals;
}

bool normal_agrees(const Vec3& triangle_normal, const Vec3& vertex_normal)
{
    if (vertex_normal.length() < 0.5)
        return true;
    return std::abs(CGLA::dot(triangle_normal, vertex_normal)) >= kNormalAgreement;
}

// Neighbor id -> (neighbor id, triangle area). One undirected chord is stored twice.
using Adj = std::map<int, std::vector<std::pair<int, double>>>;

struct LocalFan {
    bool closed = false;
    std::vector<int> order;
};

struct FanKey {
    int edges = 0;
    double area = 0;
    std::vector<int> canon;
};

bool fan_better(const FanKey& candidate, const FanKey& current)
{
    if (candidate.edges != current.edges)
        return candidate.edges > current.edges;
    if (candidate.area < current.area)
        return true;
    if (current.area < candidate.area)
        return false;
    return candidate.canon < current.canon;
}

void add_chord(Adj& adj, int a, int b, double area)
{
    if (a == b)
        return;
    for (const auto& edge : adj[a])
        if (edge.first == b)
            return;
    adj[a].push_back({b, area});
    adj[b].push_back({a, area});
}

double chord_area(const Adj& adj, int a, int b)
{
    const auto it = adj.find(a);
    if (it == adj.end())
        return 0;
    for (const auto& edge : it->second)
        if (edge.first == b)
            return edge.second;
    return 0;
}

double sequence_area(const Adj& adj, const std::vector<int>& order, bool closed)
{
    if (order.size() < 2)
        return 0;
    double area = 0;
    const size_t steps = closed ? order.size() : order.size() - 1;
    for (size_t i = 0; i < steps; ++i) {
        const int a = order[i];
        const int b = order[(i + 1) % order.size()];
        area += chord_area(adj, a, b);
    }
    return area;
}

std::vector<int> canonical_cycle(const std::vector<int>& cycle)
{
    const auto at = std::min_element(cycle.begin(), cycle.end());
    std::vector<int> forward;
    forward.reserve(cycle.size());
    for (size_t i = 0; i < cycle.size(); ++i)
        forward.push_back(cycle[(static_cast<size_t>(at - cycle.begin()) + i) % cycle.size()]);
    std::vector<int> reverse;
    reverse.reserve(forward.size());
    reverse.push_back(forward[0]);
    for (size_t i = 1; i < forward.size(); ++i)
        reverse.push_back(forward[forward.size() - i]);
    if (reverse < forward)
        return reverse;
    return forward;
}

std::vector<int> canonical_path(std::vector<int> path)
{
    std::vector<int> reverse(path.rbegin(), path.rend());
    if (reverse < path)
        return reverse;
    return path;
}

int max_degree(const Adj& adj)
{
    int degree = 0;
    for (const auto& entry : adj)
        degree = std::max(degree, static_cast<int>(entry.second.size()));
    return degree;
}

int chord_count(const Adj& adj)
{
    int count = 0;
    for (const auto& entry : adj)
        count += static_cast<int>(entry.second.size());
    return count / 2;
}

void remove_chord(Adj& adj, int a, int b)
{
    auto erase_end = [&](int from, int to) {
        const auto it = adj.find(from);
        if (it == adj.end())
            return;
        auto& edges = it->second;
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                   [&](const std::pair<int, double>& edge) { return edge.first == to; }),
                    edges.end());
        if (edges.empty())
            adj.erase(it);
    };
    erase_end(a, b);
    erase_end(b, a);
}

// Drop the largest chord at a branch until every vertex has degree at most 2.
void peel_branches(Adj& adj)
{
    const int limit = chord_count(adj) + 1;
    for (int step = 0; step < limit && max_degree(adj) > 2; ++step) {
        int best_u = -1;
        int best_v = -1;
        double best_area = -1;
        for (const auto& entry : adj) {
            if (static_cast<int>(entry.second.size()) <= 2)
                continue;
            for (const auto& edge : entry.second) {
                const int u = entry.first;
                const int v = edge.first;
                if (u > v)
                    continue;
                bool better = best_u < 0 || edge.second > best_area;
                if (!better && edge.second == best_area) {
                    const int max_id = std::max(u, v);
                    const int best_max = std::max(best_u, best_v);
                    const int min_id = std::min(u, v);
                    const int best_min = std::min(best_u, best_v);
                    better = max_id > best_max || (max_id == best_max && min_id > best_min);
                }
                if (better) {
                    best_u = u;
                    best_v = v;
                    best_area = edge.second;
                }
            }
        }
        if (best_u < 0)
            break;
        remove_chord(adj, best_u, best_v);
    }
}

std::pair<int, int> chord_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return {a, b};
}

bool order_has_chord(const std::vector<int>& order, bool closed, int a, int b)
{
    if (order.size() < 2)
        return false;
    const size_t steps = closed ? order.size() : order.size() - 1;
    for (size_t i = 0; i < steps; ++i) {
        const int u = order[i];
        const int v = order[(i + 1) % order.size()];
        if ((u == a && v == b) || (u == b && v == a))
            return true;
    }
    return false;
}

bool covers_required(const std::vector<int>& order, bool closed,
                     const std::set<std::pair<int, int>>& required)
{
    for (const auto& edge : required)
        if (!order_has_chord(order, closed, edge.first, edge.second))
            return false;
    return true;
}

bool hits_required(const std::vector<int>& order, bool closed,
                   const std::set<std::pair<int, int>>& required)
{
    for (const auto& edge : required)
        if (order_has_chord(order, closed, edge.first, edge.second))
            return true;
    return false;
}

struct CycleSearch {
    const Adj* adj = nullptr;
    const std::set<std::pair<int, int>>* required = nullptr;
    int budget = 0;
    bool exceeded = false;
    bool found = false;
    FanKey best;
};

void consider_cycle(CycleSearch& search, const std::vector<int>& cycle)
{
    if (cycle.size() < 3)
        return;
    if (search.required != nullptr && !covers_required(cycle, true, *search.required))
        return;
    FanKey key;
    key.edges = static_cast<int>(cycle.size());
    key.area = sequence_area(*search.adj, cycle, true);
    key.canon = canonical_cycle(cycle);
    if (!search.found || fan_better(key, search.best)) {
        search.found = true;
        search.best = std::move(key);
    }
}

void search_cycles_from(CycleSearch& search, int start, std::vector<int>& path, std::map<int, char>& on_path)
{
    if (search.budget-- <= 0) {
        search.exceeded = true;
        return;
    }
    const int at = path.back();
    for (const auto& edge : search.adj->at(at)) {
        const int next = edge.first;
        if (next == start && path.size() >= 3) {
            consider_cycle(search, path);
            continue;
        }
        if (next <= start || on_path[next])
            continue;
        path.push_back(next);
        on_path[next] = 1;
        search_cycles_from(search, start, path, on_path);
        if (search.exceeded)
            return;
        on_path[next] = 0;
        path.pop_back();
    }
}

std::vector<int> component_nodes(const Adj& adj, int seed, std::map<int, char>& seen)
{
    std::vector<int> nodes;
    std::vector<int> stack{seed};
    seen[seed] = 1;
    while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        nodes.push_back(node);
        for (const auto& edge : adj.at(node)) {
            if (!seen[edge.first]) {
                seen[edge.first] = 1;
                stack.push_back(edge.first);
            }
        }
    }
    return nodes;
}

std::vector<int> walk_chain(const Adj& adj, int start, int count)
{
    std::vector<int> order;
    int previous = -1;
    int current = start;
    for (int step = 0; step < count; ++step) {
        order.push_back(current);
        int next = -1;
        for (const auto& edge : adj.at(current)) {
            if (edge.first != previous) {
                next = edge.first;
                break;
            }
        }
        if (next < 0)
            break;
        previous = current;
        current = next;
        if (current == start)
            break;
    }
    return order;
}

std::vector<int> tree_path(const Adj& adj, int from, int to)
{
    std::map<int, int> parent;
    parent[from] = -1;
    std::vector<int> stack{from};
    while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        for (const auto& edge : adj.at(node)) {
            if (parent.find(edge.first) != parent.end())
                continue;
            parent[edge.first] = node;
            stack.push_back(edge.first);
        }
    }
    std::vector<int> path;
    const auto destination = parent.find(to);
    if (destination == parent.end())
        return path;
    for (int current = to; current >= 0; current = parent[current])
        path.push_back(current);
    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<LocalFan> longest_paths(const Adj& adj)
{
    std::vector<LocalFan> fans;
    std::map<int, char> seen;
    for (const auto& entry : adj) {
        if (seen[entry.first])
            continue;
        const std::vector<int> nodes = component_nodes(adj, entry.first, seen);
        std::vector<int> leaves;
        for (int node : nodes)
            if (adj.at(node).size() <= 1)
                leaves.push_back(node);
        if (leaves.size() < 2)
            continue;
        FanKey best;
        bool any = false;
        for (size_t i = 0; i < leaves.size(); ++i) {
            for (size_t j = i + 1; j < leaves.size(); ++j) {
                const std::vector<int> path = tree_path(adj, leaves[i], leaves[j]);
                if (path.size() < 2)
                    continue;
                FanKey key;
                key.edges = static_cast<int>(path.size()) - 1;
                key.area = sequence_area(adj, path, false);
                key.canon = canonical_path(path);
                if (!any || fan_better(key, best)) {
                    any = true;
                    best = std::move(key);
                }
            }
        }
        if (any)
            fans.push_back(LocalFan{false, std::move(best.canon)});
    }
    return fans;
}

struct Chains {
    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<int>> paths;
};

Chains collect_chains(const Adj& adj)
{
    Chains chains;
    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<int>> paths;
    std::map<int, char> seen;
    for (const auto& entry : adj) {
        if (seen[entry.first])
            continue;
        const std::vector<int> nodes = component_nodes(adj, entry.first, seen);
        if (nodes.size() < 2)
            continue;
        std::vector<int> ends;
        int degree = 0;
        for (int node : nodes) {
            const int node_degree = static_cast<int>(adj.at(node).size());
            degree = std::max(degree, node_degree);
            if (node_degree <= 1)
                ends.push_back(node);
        }
        if (degree > 2)
            continue;
        if (ends.empty() && nodes.size() >= 3) {
            const int start = *std::min_element(nodes.begin(), nodes.end());
            std::vector<int> cycle = walk_chain(adj, start, static_cast<int>(nodes.size()));
            int previous = cycle.size() >= 2 ? cycle[cycle.size() - 2] : -1;
            int next = -1;
            if (!cycle.empty()) {
                for (const auto& edge : adj.at(cycle.back())) {
                    if (edge.first != previous) {
                        next = edge.first;
                        break;
                    }
                }
            }
            if (cycle.size() == nodes.size() && next == start)
                cycles.push_back(std::move(cycle));
            continue;
        }
        if (ends.empty())
            continue;
        const int start = *std::min_element(ends.begin(), ends.end());
        std::vector<int> path = walk_chain(adj, start, static_cast<int>(nodes.size()));
        if (path.size() >= 2)
            paths.push_back(std::move(path));
    }
    chains.cycles = std::move(cycles);
    chains.paths = std::move(paths);
    return chains;
}

std::vector<LocalFan> fans_from_chains(const Adj& adj)
{
    const Chains chains = collect_chains(adj);
    if (!chains.cycles.empty()) {
        FanKey best;
        bool any = false;
        for (const auto& cycle : chains.cycles) {
            FanKey key;
            key.edges = static_cast<int>(cycle.size());
            key.area = sequence_area(adj, cycle, true);
            key.canon = canonical_cycle(cycle);
            if (!any || fan_better(key, best)) {
                any = true;
                best = std::move(key);
            }
        }
        if (any)
            return {LocalFan{true, std::move(best.canon)}};
    }

    std::vector<LocalFan> fans;
    for (const auto& path : chains.paths) {
        std::vector<int> canon = canonical_path(path);
        if (canon.size() >= 2)
            fans.push_back(LocalFan{false, std::move(canon)});
    }
    return fans;
}

bool chord_present(const Adj& adj, int a, int b)
{
    const auto it = adj.find(a);
    if (it == adj.end())
        return false;
    for (const auto& edge : it->second)
        if (edge.first == b)
            return true;
    return false;
}

// Keep committed chords. A cycle that contains every one of them wins.
// Otherwise each remaining chain that touches one is its own fan.
std::vector<LocalFan> fans_covering(const Adj& adj, const std::set<std::pair<int, int>>& required)
{
    if (max_degree(adj) > 2) {
        std::vector<LocalFan> fans;
        for (const auto& edge : required) {
            if (!chord_present(adj, edge.first, edge.second))
                continue;
            fans.push_back(LocalFan{false, canonical_path({edge.first, edge.second})});
        }
        return fans;
    }

    const Chains chains = collect_chains(adj);
    FanKey best;
    bool covered = false;
    for (const auto& cycle : chains.cycles) {
        if (!covers_required(cycle, true, required))
            continue;
        FanKey key;
        key.edges = static_cast<int>(cycle.size());
        key.area = sequence_area(adj, cycle, true);
        key.canon = canonical_cycle(cycle);
        if (!covered || fan_better(key, best)) {
            covered = true;
            best = std::move(key);
        }
    }
    if (covered)
        return {LocalFan{true, std::move(best.canon)}};

    std::vector<LocalFan> fans;
    for (const auto& cycle : chains.cycles) {
        if (!hits_required(cycle, true, required))
            continue;
        fans.push_back(LocalFan{true, canonical_cycle(cycle)});
    }
    for (const auto& path : chains.paths) {
        if (!hits_required(path, false, required))
            continue;
        std::vector<int> canon = canonical_path(path);
        if (canon.size() >= 2)
            fans.push_back(LocalFan{false, std::move(canon)});
    }
    return fans;
}

void peel_non_required(Adj& adj, const std::set<std::pair<int, int>>& required)
{
    const int limit = chord_count(adj) + 1;
    for (int step = 0; step < limit && max_degree(adj) > 2; ++step) {
        int best_u = -1;
        int best_v = -1;
        double best_area = -1;
        for (const auto& entry : adj) {
            if (static_cast<int>(entry.second.size()) <= 2)
                continue;
            for (const auto& edge : entry.second) {
                const int u = entry.first;
                const int v = edge.first;
                if (u > v || required.count(chord_key(u, v)) != 0)
                    continue;
                bool better = best_u < 0 || edge.second > best_area;
                if (!better && edge.second == best_area) {
                    const int max_id = std::max(u, v);
                    const int best_max = std::max(best_u, best_v);
                    const int min_id = std::min(u, v);
                    const int best_min = std::min(best_u, best_v);
                    better = max_id > best_max || (max_id == best_max && min_id > best_min);
                }
                if (better) {
                    best_u = u;
                    best_v = v;
                    best_area = edge.second;
                }
            }
        }
        if (best_u < 0)
            break;
        remove_chord(adj, best_u, best_v);
    }
}

// A few dozen chords is the enumeration limit. Past that, peeling is linear.
constexpr int kCycleSearchEdges = 48;
constexpr int kCycleSearchBudget = 200000;

std::vector<LocalFan> search_longest_cycle(const Adj& adj, const std::set<std::pair<int, int>>* required,
                                           bool& exceeded)
{
    exceeded = false;
    CycleSearch search;
    search.adj = &adj;
    search.required = required;
    search.budget = kCycleSearchBudget;
    for (const auto& entry : adj) {
        std::vector<int> path{entry.first};
        std::map<int, char> on_path;
        on_path[entry.first] = 1;
        search_cycles_from(search, entry.first, path, on_path);
        if (search.exceeded)
            break;
    }
    exceeded = search.exceeded;
    if (search.exceeded || !search.found)
        return {};
    return {LocalFan{true, std::move(search.best.canon)}};
}

std::vector<LocalFan> choose_fans(Adj adj, const std::set<std::pair<int, int>>& required)
{
    if (!required.empty()) {
        if (chord_count(adj) > kCycleSearchEdges) {
            peel_non_required(adj, required);
            return fans_covering(adj, required);
        }
        bool exceeded = false;
        std::vector<LocalFan> cycle = search_longest_cycle(adj, &required, exceeded);
        if (!exceeded && !cycle.empty())
            return cycle;
        peel_non_required(adj, required);
        return fans_covering(adj, required);
    }

    if (chord_count(adj) > kCycleSearchEdges) {
        peel_branches(adj);
    } else if (max_degree(adj) > 2) {
        bool exceeded = false;
        std::vector<LocalFan> cycle = search_longest_cycle(adj, nullptr, exceeded);
        if (exceeded)
            peel_branches(adj);
        else if (!cycle.empty())
            return cycle;
    }

    if (max_degree(adj) <= 2)
        return fans_from_chains(adj);
    return longest_paths(adj);
}

// Both endpoints fit in 32 bits. reconstruct rejects a larger cloud.
uint64_t undirected_id(size_t a, size_t b)
{
    if (a > b)
        std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
}

uint64_t directed_id(size_t from, size_t to)
{
    return (static_cast<uint64_t>(from) << 32) | static_cast<uint32_t>(to);
}

// A triangle this close to the plane of the edge and its normal is neither side.
constexpr double kSideEpsilon = 1e-8;

bool finite_vector(const Vec3& value)
{
    return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

// One Delaunay triangle on each side of every tunnel-weighted forest edge.
// Left of directed edge a→b is the side where (b-a) × (c-a) agrees with the
// edge normal. The open side keeps the triangle of least
// perimeter * area / projected area.
// A triangle already chosen for an earlier forest edge, with a second edge
// on this forest edge, fills that side and is not replaced. Every edge of a
// chosen triangle records the side that triangle occupies. A new candidate
// on an occupied side of any of its edges is skipped, and the next score is
// tried. A candidate that would be the third triangle on any of its edges is
// skipped the same way. So is one whose spoke at a vertex falls strictly
// between the edges of a triangle already incident there, measured in the
// plane of that vertex normal. A spoke that is an edge of that triangle is
// not tested against it. The stored winding puts
// the edge normal on a newly chosen triangle, so the right triangle uses the
// opposite direction along the edge.
void add_forest_triangles(const std::vector<Vec3>& points, const LRsROpts& opts, Manifold& result)
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

    std::unordered_map<uint64_t, char> forest_edge;
    for (Geometry::AMGraph::NodeID node : forest.node_ids()) {
        for (const auto& entry : forest.edges(node)) {
            if (node < entry.first)
                forest_edge.emplace(undirected_id(node, entry.first), 1);
        }
    }

    std::unordered_map<TriKey, char, TriKeyHash> have;
    std::vector<std::array<size_t, 3>> faces;
    std::unordered_map<uint64_t, std::vector<size_t>> faces_on;
    // Bit 1 is the left side of an edge, bit 2 the right side.
    std::unordered_map<uint64_t, unsigned char> side_taken;
    // The two neighbor vertices of each triangle already chosen at a vertex.
    std::vector<std::vector<std::pair<size_t, size_t>>> vertex_tris(points.size());
    auto edge_side = [&](size_t u, size_t v, size_t corner) {
        if (u > v)
            std::swap(u, v);
        if (corner == u || corner == v)
            return 0;
        Vec3 normal_u = normals[u];
        Vec3 normal_v = normals[v];
        if (!finite_vector(normal_u) || !finite_vector(normal_v))
            return 0;
        if (CGLA::dot(normal_u, normal_v) < 0.0)
            normal_v = -normal_v;
        const Vec3 edge_normal = normal_u + normal_v;
        const double normal_length = edge_normal.length();
        if (!(normal_length > 1e-12) || !std::isfinite(normal_length))
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
    };
    auto mark_side = [&](size_t u, size_t v, int side) {
        if (side == 0)
            return;
        side_taken[undirected_id(u, v)] |= side > 0 ? 1 : 2;
    };
    auto side_occupied = [&](size_t u, size_t v, int side) {
        if (side == 0)
            return false;
        const auto found = side_taken.find(undirected_id(u, v));
        if (found == side_taken.end())
            return false;
        const unsigned char bit = side > 0 ? 1 : 2;
        return (found->second & bit) != 0;
    };
    // The candidate occupies one side of edge (u, v). Reject it when that side
    // already holds a triangle. Forest edges and other edges use this same test.
    auto side_blocked = [&](size_t u, size_t v, size_t other) {
        return side_occupied(u, v, edge_side(u, v, other));
    };
    auto incident_faces = [&](size_t u, size_t v) {
        const auto found = faces_on.find(undirected_id(u, v));
        if (found == faces_on.end())
            return 0;
        return static_cast<int>(found->second.size());
    };
    auto emit = [&](size_t i, size_t j, size_t k) {
        if (!have.emplace(tri_key(i, j, k), 1).second)
            return;
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
    };
    // Project the edges onto the plane perpendicular to the vertex normal.
    // The spoke is between the other two when that projection is a positive
    // combination of theirs. Flipping the normal leaves the test unchanged.
    // An angle that collapses in the tangent plane does not block, and the
    // boundary of the angle is not inside.
    auto lies_between = [&](size_t vertex, size_t end_a, size_t end_b, size_t spoke) {
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
    };
    auto vertex_blocked = [&](size_t vertex, size_t spoke_u, size_t spoke_v) {
        for (const auto& triangle : vertex_tris[vertex]) {
            const size_t end_a = triangle.first;
            const size_t end_b = triangle.second;
            if (spoke_u != end_a && spoke_u != end_b && lies_between(vertex, end_a, end_b, spoke_u))
                return true;
            if (spoke_v != end_a && spoke_v != end_b && lies_between(vertex, end_a, end_b, spoke_v))
                return true;
        }
        return false;
    };
    auto forest_edges_of = [&](const std::array<size_t, 3>& face) {
        int count = 0;
        if (forest_edge.count(undirected_id(face[0], face[1])) != 0)
            ++count;
        if (forest_edge.count(undirected_id(face[1], face[2])) != 0)
            ++count;
        if (forest_edge.count(undirected_id(face[2], face[0])) != 0)
            ++count;
        return count;
    };

    for (Geometry::AMGraph::NodeID a : forest.node_ids()) {
        for (const auto& entry : forest.edges(a)) {
            const Geometry::AMGraph::NodeID b = entry.first;
            if (a >= b)
                continue;
            Vec3 normal_a = normals[a];
            Vec3 normal_b = normals[b];
            if (!finite_vector(normal_a) || !finite_vector(normal_b))
                continue;
            // The higher index is the end that flips when the normals disagree.
            if (CGLA::dot(normal_a, normal_b) < 0.0)
                normal_b = -normal_b;
            const Vec3 edge_normal = normal_a + normal_b;
            const double normal_length = edge_normal.length();
            if (!(normal_length > 1e-12) || !std::isfinite(normal_length))
                continue;

            const Vec3& point_a = points[a];
            const Vec3& point_b = points[b];
            const Vec3 along = point_b - point_a;
            // Earlier forest edges are the ones already visited. A chosen
            // triangle that also contains this edge occupies its side.
            bool left_taken = false;
            bool right_taken = false;
            if (const auto found = faces_on.find(undirected_id(a, b)); found != faces_on.end()) {
                for (size_t face_index : found->second) {
                    const std::array<size_t, 3>& face = faces[face_index];
                    if (forest_edges_of(face) < 2)
                        continue;
                    size_t corner = face[0];
                    if (corner == a || corner == b)
                        corner = face[1];
                    if (corner == a || corner == b)
                        corner = face[2];
                    const double area = triangle_area(point_a, point_b, points[corner]);
                    if (!(area > 0.0) || !std::isfinite(area))
                        continue;
                    const double side = CGLA::dot(
                        CGLA::cross(along, points[corner] - point_a), edge_normal);
                    if (!std::isfinite(side))
                        continue;
                    const double cosine = side / (2.0 * area * normal_length);
                    if (cosine > kSideEpsilon)
                        left_taken = true;
                    else if (cosine < -kSideEpsilon)
                        right_taken = true;
                }
            }
            if (left_taken)
                mark_side(a, b, 1);
            if (right_taken)
                mark_side(a, b, -1);
            if (left_taken && right_taken)
                continue;

            int left = -1;
            int right = -1;
            double left_score = 0.0;
            double right_score = 0.0;
            for (Geometry::AMGraph::NodeID corner : graph.shared_neighbors(a, b)) {
                if (corner == a || corner == b)
                    continue;
                const Vec3& point_c = points[corner];
                const double area = triangle_area(point_a, point_b, point_c);
                if (!(area > 0.0) || !std::isfinite(area))
                    continue;
                const double side = CGLA::dot(CGLA::cross(along, point_c - point_a), edge_normal);
                if (!std::isfinite(side))
                    continue;
                const double cosine = side / (2.0 * area * normal_length);
                if (!(std::abs(cosine) > kSideEpsilon))
                    continue;
                if ((cosine > 0.0 && left_taken) || (cosine < 0.0 && right_taken))
                    continue;
                // An occupied side of any edge drops this candidate. The next
                // score on this side of the forest edge is tried.
                if (side_blocked(a, b, corner) || side_blocked(a, corner, b)
                    || side_blocked(b, corner, a))
                    continue;
                // Two faces already meet on one of these edges. Another face there
                // would be the third, so the next score is tried instead.
                if (incident_faces(a, b) >= 2 || incident_faces(a, corner) >= 2
                    || incident_faces(b, corner) >= 2)
                    continue;
                // A non-shared spoke whose tangent, in the vertex normal's plane,
                // falls inside an incident triangle's angle drops the candidate.
                // The next score is tried instead.
                if (vertex_blocked(a, b, corner) || vertex_blocked(b, a, corner)
                    || vertex_blocked(corner, a, b))
                    continue;
                const double score = projected_perimeter_score(point_a, point_b, point_c, edge_normal);
                if (!std::isfinite(score))
                    continue;
                int& chosen = cosine > 0.0 ? left : right;
                double& chosen_score = cosine > 0.0 ? left_score : right_score;
                const int index = static_cast<int>(corner);
                if (chosen >= 0 && (chosen_score < score || (chosen_score == score && chosen <= index)))
                    continue;
                chosen = index;
                chosen_score = score;
            }
            // Both new triangles use this forest edge. One existing face leaves
            // room for only one of them.
            if (left >= 0 && right >= 0 && incident_faces(a, b) >= 1) {
                const bool keep_left = left_score < right_score
                    || (left_score == right_score && left <= right);
                if (keep_left)
                    right = -1;
                else
                    left = -1;
            }
            if (left >= 0) {
                emit(a, b, static_cast<size_t>(left));
                mark_side(a, b, 1);
            }
            if (right >= 0) {
                emit(a, static_cast<size_t>(right), b);
                mark_side(a, b, -1);
            }
        }
    }

    for (const auto& face : faces) {
        const std::array<Vec3, 3> corners{
            points[face[0]], points[face[1]], points[face[2]]};
        result.add_face(corners);
    }
    stitch_mesh(result, 1e-8);
    result.cleanup();
}

} // namespace

void reconstruct(const std::vector<Vec3>& points, const LRsROpts& opts, Manifold& result)
{
    result.clear();
    if (points.empty() || points.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        return;

    // Longest-cycle fans at each vertex stay compiled below. The forest walk is active.
    constexpr bool kVertexFans = false;
    if (!kVertexFans) {
        add_forest_triangles(points, opts, result);
        return;
    }

    const int n = static_cast<int>(points.size());

    const double spacing = median_neighbor_spacing(points);
    const double radius = (opts.radius_factor > 0.0) ? opts.radius_factor * spacing : 0.0;
    const std::vector<Vec3> normals = local_normals(points);
    const std::vector<Geometry::DelaunayStar> stars = Geometry::delaunay_stars(points, radius);
    if (static_cast<int>(stars.size()) != n)
        return;

    std::unordered_map<uint64_t, char> graph_edges;
    auto edge_id = [](int a, int b) {
        if (a > b)
            std::swap(a, b);
        return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32)
            | static_cast<uint32_t>(b);
    };
    auto is_edge = [&](int a, int b) {
        if (a < 0 || b < 0 || a >= n || b >= n || a == b)
            return false;
        return graph_edges.find(edge_id(a, b)) != graph_edges.end();
    };
    for (int i = 0; i < n; ++i) {
        for (int neighbor : stars[static_cast<size_t>(i)].neighbors) {
            if (neighbor < 0 || neighbor >= n || neighbor == i)
                continue;
            graph_edges.emplace(edge_id(i, neighbor), 1);
        }
    }

    // At vertex v the chord (a, b) is already committed as face (v, a, b).
    std::vector<std::vector<std::pair<int, int>>> committed_chords(static_cast<size_t>(n));
    std::unordered_map<TriKey, char, TriKeyHash> have;
    std::unordered_map<uint64_t, int> edge_faces;
    std::unordered_map<uint64_t, char> directed_used;
    std::vector<std::array<size_t, 3>> faces;
    std::vector<char> placed(static_cast<size_t>(n), 0);
    std::set<int> pending;
    // One sheet grows from the lowest index that still has no triangle.
    // Committing a face queues its corners so the next fan matches that winding.
    auto enqueue = [&](int corner) {
        if (corner < 0 || corner >= n || placed[static_cast<size_t>(corner)] != 0)
            return;
        pending.insert(corner);
    };

    auto occupied = [&](const std::array<size_t, 3>& face) {
        for (int edge = 0; edge < 3; ++edge) {
            const size_t from = face[static_cast<size_t>(edge)];
            const size_t to = face[static_cast<size_t>((edge + 1) % 3)];
            if (directed_used.count(directed_id(from, to)) != 0)
                return true;
            if (edge_faces[undirected_id(from, to)] >= 2)
                return true;
        }
        return false;
    };
    auto commit_face = [&](const std::array<size_t, 3>& face) {
        const TriKey key = tri_key(face[0], face[1], face[2]);
        if (have.count(key) != 0 || occupied(face))
            return;
        have.emplace(key, 1);
        faces.push_back(face);
        for (int edge = 0; edge < 3; ++edge) {
            const size_t from = face[static_cast<size_t>(edge)];
            const size_t to = face[static_cast<size_t>((edge + 1) % 3)];
            directed_used.emplace(directed_id(from, to), 1);
            edge_faces[undirected_id(from, to)] += 1;
        }
        committed_chords[face[0]].push_back({static_cast<int>(face[1]), static_cast<int>(face[2])});
        committed_chords[face[1]].push_back({static_cast<int>(face[2]), static_cast<int>(face[0])});
        committed_chords[face[2]].push_back({static_cast<int>(face[0]), static_cast<int>(face[1])});
        enqueue(static_cast<int>(face[0]));
        enqueue(static_cast<int>(face[1]));
        enqueue(static_cast<int>(face[2]));
    };

    int cursor = 0;
    while (true) {
        if (pending.empty()) {
            while (cursor < n && placed[static_cast<size_t>(cursor)] != 0)
                ++cursor;
            if (cursor >= n)
                break;
            pending.insert(cursor);
        }
        const int vertex = *pending.begin();
        pending.erase(pending.begin());
        if (placed[static_cast<size_t>(vertex)] != 0)
            continue;
        placed[static_cast<size_t>(vertex)] = 1;

        Adj adj;
        for (const Geometry::DelaunayLink& link : stars[static_cast<size_t>(vertex)].links) {
            const int k = static_cast<int>(link.neighbors.size());
            if (k < 2)
                continue;
            const int steps = link.closed ? k : k - 1;
            for (int t = 0; t < steps; ++t) {
                const int a = link.neighbors[static_cast<size_t>(t)];
                const int b = link.neighbors[static_cast<size_t>(link.closed ? (t + 1) % k : t + 1)];
                if (!is_edge(vertex, a) || !is_edge(vertex, b) || !is_edge(a, b))
                    continue;
                const Vec3& origin = points[static_cast<size_t>(vertex)];
                const Vec3& corner_a = points[static_cast<size_t>(a)];
                const Vec3& corner_b = points[static_cast<size_t>(b)];
                const double area = triangle_area(origin, corner_a, corner_b);
                if (!(area > 0.0) || !std::isfinite(area))
                    continue;
                Vec3 triangle_normal = CGLA::cross(corner_a - origin, corner_b - origin);
                const double normal_length = triangle_normal.length();
                if (!(normal_length > 0.0) || !std::isfinite(normal_length))
                    continue;
                triangle_normal /= normal_length;
                if (!normal_agrees(triangle_normal, normals[static_cast<size_t>(vertex)])
                    || !normal_agrees(triangle_normal, normals[static_cast<size_t>(a)])
                    || !normal_agrees(triangle_normal, normals[static_cast<size_t>(b)]))
                    continue;
                add_chord(adj, a, b, area);
            }
        }
        for (auto& entry : adj)
            std::sort(entry.second.begin(), entry.second.end());

        const bool touched = !committed_chords[static_cast<size_t>(vertex)].empty();
        std::set<std::pair<int, int>> required;
        std::vector<std::pair<int, int>> wanted;
        for (const auto& chord : committed_chords[static_cast<size_t>(vertex)]) {
            if (!chord_present(adj, chord.first, chord.second))
                continue;
            required.insert(chord_key(chord.first, chord.second));
            wanted.push_back(chord);
        }
        // Committed triangles already cover this vertex. A chord that is not
        // admissible here cannot be extended, so a new fan would be disjoint.
        if (touched && required.empty())
            continue;

        std::vector<LocalFan> chosen = choose_fans(std::move(adj), required);
        for (LocalFan& fan : chosen) {
            if (fan.order.size() < 2)
                continue;
            int agree = 0;
            int disagree = 0;
            const size_t steps = fan.closed ? fan.order.size() : fan.order.size() - 1;
            for (size_t t = 0; t < steps; ++t) {
                const int a = fan.order[t];
                const int b = fan.order[(t + 1) % fan.order.size()];
                bool known = false;
                for (const auto& chord : wanted) {
                    if (chord.first == a && chord.second == b) {
                        ++agree;
                        known = true;
                    } else if (chord.first == b && chord.second == a) {
                        ++disagree;
                        known = true;
                    }
                }
                if (known)
                    continue;
                // The other face of this chord does not touch this vertex.
                // Use the free direction so the two faces sit on opposite sides.
                const size_t ua = static_cast<size_t>(a);
                const size_t ub = static_cast<size_t>(b);
                if (directed_used.count(directed_id(ua, ub)) != 0)
                    ++disagree;
                else if (directed_used.count(directed_id(ub, ua)) != 0)
                    ++agree;
            }
            if (disagree > agree)
                std::reverse(fan.order.begin(), fan.order.end());
            const size_t oriented_steps = fan.closed ? fan.order.size() : fan.order.size() - 1;
            for (size_t t = 0; t < oriented_steps; ++t) {
                const size_t a = static_cast<size_t>(fan.order[t]);
                const size_t b = static_cast<size_t>(fan.order[(t + 1) % fan.order.size()]);
                commit_face({static_cast<size_t>(vertex), a, b});
            }
        }
    }

    for (const auto& face : faces) {
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
