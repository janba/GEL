/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Delaunay3D.h"

#include <GEL/Geometry/KDTree.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Geometry {
namespace {

using CGLA::Vec3d;
using CGLA::cross;
using CGLA::dot;
using CGLA::normalize;
using CGLA::sqr_length;

struct Plane {
    Vec3d n;
    double d = 0;
    int src = -1;
};

struct Vert {
    Vec3d p;
    std::vector<int> planes;
};

struct Edge {
    int a = 0;
    int b = 0;
    int p0 = -1;
    int p1 = -1;
};

struct Cell {
    std::vector<Plane> planes;
    std::vector<Vert> verts;
    std::vector<Edge> edges;
    double half = 0;
    double tol = 0;
};

enum class Side { In, On, Out };

double side_of(const Plane& pl, const Vec3d& x)
{
    return dot(pl.n, x) + pl.d;
}

Side classify(double s, double tol)
{
    if (s > tol)
        return Side::Out;
    if (s < -tol)
        return Side::In;
    return Side::On;
}

bool finite_point(const Vec3d& x)
{
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

bool intersect_planes(const Plane& a, const Plane& b, const Plane& c, Vec3d& x)
{
    const Vec3d bxc = cross(b.n, c.n);
    const double det = dot(a.n, bxc);
    if (!(std::abs(det) > 1e-14))
        return false;
    const Vec3d cxa = cross(c.n, a.n);
    const Vec3d axb = cross(a.n, b.n);
    x = (bxc * (-a.d) + cxa * (-b.d) + axb * (-c.d)) / det;
    return finite_point(x);
}

bool interpolate_edge(const Plane& pl, const Vec3d& a, const Vec3d& b, Vec3d& x)
{
    const double s0 = side_of(pl, a);
    const double s1 = side_of(pl, b);
    const double denom = s0 - s1;
    if (!(std::abs(denom) > 0))
        return false;
    const double t = s0 / denom;
    if (t < -1e-8 || t > 1.0 + 1e-8)
        return false;
    x = a * (1.0 - t) + b * t;
    return finite_point(x);
}

bool inside_cell(const Cell& cell, const Vec3d& x)
{
    const double limit = cell.tol * 8.0;
    for (const Plane& pl : cell.planes)
        if (side_of(pl, x) > limit)
            return false;
    return true;
}

void plane_axes(const Vec3d& n, Vec3d& u, Vec3d& v)
{
    const Vec3d axis = std::abs(n[2]) < 0.9 ? Vec3d(0, 0, 1) : Vec3d(0, 1, 0);
    u = normalize(cross(n, axis));
    v = cross(n, u);
}

bool has_edge(const std::vector<Edge>& edges, int a, int b)
{
    for (const Edge& e : edges)
        if ((e.a == a && e.b == b) || (e.a == b && e.b == a))
            return true;
    return false;
}

int shared_plane(const Vert& a, const Vert& b, int fresh, const std::vector<Plane>& planes, double tol)
{
    for (int pa : a.planes)
        for (int pb : b.planes)
            if (pa == pb && pa != fresh)
                return pa;
    const double limit = tol * 32.0;
    for (int pi = 0; pi < static_cast<int>(planes.size()); ++pi) {
        if (pi == fresh)
            continue;
        if (std::abs(side_of(planes[pi], a.p)) <= limit &&
            std::abs(side_of(planes[pi], b.p)) <= limit)
            return pi;
    }
    return -1;
}

double polygon_cross2(const std::vector<Vert>& verts, const std::vector<int>& ids)
{
    if (ids.size() < 3)
        return 0;
    const Vec3d origin = verts[ids[0]].p;
    double best = 0;
    for (size_t i = 1; i < ids.size(); ++i) {
        for (size_t j = i + 1; j < ids.size(); ++j) {
            const double c2 = sqr_length(cross(verts[ids[i]].p - origin, verts[ids[j]].p - origin));
            if (c2 > best)
                best = c2;
        }
    }
    return best;
}

/// Cube [-half, half]^3 in site-local coordinates. Inside is n·x + d <= 0.
void init_cube(Cell& cell, double half)
{
    cell.half = half;
    cell.tol = std::max(1e-15, 1e-10 * half);
    cell.planes = {
        {Vec3d(-1, 0, 0), -half, -1},
        {Vec3d( 1, 0, 0), -half, -1},
        {Vec3d(0, -1, 0), -half, -1},
        {Vec3d(0,  1, 0), -half, -1},
        {Vec3d(0, 0, -1), -half, -1},
        {Vec3d(0, 0,  1), -half, -1},
    };
    const double h = half;
    const Vec3d corners[8] = {
        Vec3d(-h, -h, -h), Vec3d( h, -h, -h), Vec3d(-h,  h, -h), Vec3d( h,  h, -h),
        Vec3d(-h, -h,  h), Vec3d( h, -h,  h), Vec3d(-h,  h,  h), Vec3d( h,  h,  h),
    };
    const int corner_planes[8][3] = {
        {0, 2, 4}, {1, 2, 4}, {0, 3, 4}, {1, 3, 4},
        {0, 2, 5}, {1, 2, 5}, {0, 3, 5}, {1, 3, 5},
    };
    cell.verts.resize(8);
    for (int i = 0; i < 8; ++i) {
        cell.verts[i].p = corners[i];
        cell.verts[i].planes = {corner_planes[i][0], corner_planes[i][1], corner_planes[i][2]};
    }
    const int cube_edges[12][4] = {
        {0, 1, 2, 4}, {2, 3, 3, 4}, {4, 5, 2, 5}, {6, 7, 3, 5},
        {0, 2, 0, 4}, {1, 3, 1, 4}, {4, 6, 0, 5}, {5, 7, 1, 5},
        {0, 4, 0, 2}, {1, 5, 1, 2}, {2, 6, 0, 3}, {3, 7, 1, 3},
    };
    cell.edges.resize(12);
    for (int i = 0; i < 12; ++i) {
        cell.edges[i].a = cube_edges[i][0];
        cell.edges[i].b = cube_edges[i][1];
        cell.edges[i].p0 = cube_edges[i][2];
        cell.edges[i].p1 = cube_edges[i][3];
    }
}

bool misses_box(const Cell& cell, const Plane& pl)
{
    Vec3d lo = cell.verts[0].p;
    Vec3d hi = cell.verts[0].p;
    for (const Vert& v : cell.verts) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], v.p[k]);
            hi[k] = std::max(hi[k], v.p[k]);
        }
    }
    double max_side = pl.d;
    for (int k = 0; k < 3; ++k)
        max_side += pl.n[k] * (pl.n[k] >= 0 ? hi[k] : lo[k]);
    return max_side < -cell.tol;
}

bool crossing_point(const Cell& cell, const Plane& pl, const Edge& e, Vec3d& x)
{
    Vec3d linear;
    const bool have_linear = interpolate_edge(pl, cell.verts[e.a].p, cell.verts[e.b].p, linear);
    if (e.p0 >= 0 && e.p1 >= 0) {
        Vec3d exact;
        if (intersect_planes(cell.planes[e.p0], cell.planes[e.p1], pl, exact) && inside_cell(cell, exact)) {
            if (!have_linear) {
                x = exact;
                return true;
            }
            const double slack = std::max(cell.tol * 100.0, 1e-8 * cell.half);
            if (sqr_length(exact - linear) <= slack * slack) {
                x = exact;
                return true;
            }
        }
    }
    if (!have_linear)
        return false;
    x = linear;
    return true;
}

/**
 * Clip the convex cell by one bisector. A crossed edge contributes the two
 * planes stored on that edge; the new vertex is their intersection with `pl`.
 * Returns false when the plane misses the cell.
 */
bool clip_by_plane(Cell& cell, const Plane& pl)
{
    if (cell.verts.empty() || misses_box(cell, pl))
        return false;

    const int nv = static_cast<int>(cell.verts.size());
    std::vector<Side> cls(nv);
    bool any_out = false;
    for (int i = 0; i < nv; ++i) {
        cls[i] = classify(side_of(pl, cell.verts[i].p), cell.tol);
        if (cls[i] == Side::Out)
            any_out = true;
    }
    if (!any_out)
        return false;

    const int fresh = static_cast<int>(cell.planes.size());
    std::vector<Vert> verts;
    verts.reserve(cell.verts.size() + 8);
    std::vector<int> remap(nv, -1);
    for (int i = 0; i < nv; ++i) {
        if (cls[i] == Side::Out)
            continue;
        remap[i] = static_cast<int>(verts.size());
        verts.push_back(cell.verts[i]);
        if (cls[i] == Side::On)
            verts.back().planes.push_back(fresh);
    }

    std::vector<Edge> edges;
    edges.reserve(cell.edges.size() + 8);
    std::vector<int> cut;
    cut.reserve(16);
    std::vector<char> in_cut(nv, 0);

    auto mark_old = [&](int old_i) {
        if (!in_cut[old_i]) {
            in_cut[old_i] = 1;
            cut.push_back(remap[old_i]);
        }
    };

    for (const Edge& e : cell.edges) {
        const Side sa = cls[e.a];
        const Side sb = cls[e.b];
        if (sa == Side::Out && sb == Side::Out)
            continue;
        if (sa != Side::Out && sb != Side::Out) {
            Edge kept = e;
            kept.a = remap[e.a];
            kept.b = remap[e.b];
            edges.push_back(kept);
            continue;
        }

        const int in_old = sa == Side::Out ? e.b : e.a;
        if (cls[in_old] == Side::On) {
            mark_old(in_old);
            continue;
        }

        Vec3d x;
        if (!crossing_point(cell, pl, e, x))
            return false;

        int snap = -1;
        const double snap2 = cell.tol * cell.tol;
        for (int i = 0; i < static_cast<int>(verts.size()); ++i) {
            if (sqr_length(verts[i].p - x) <= snap2) {
                snap = i;
                break;
            }
        }
        int created = snap;
        if (created < 0) {
            Vert nvx;
            nvx.p = x;
            if (e.p0 >= 0)
                nvx.planes.push_back(e.p0);
            if (e.p1 >= 0)
                nvx.planes.push_back(e.p1);
            nvx.planes.push_back(fresh);
            created = static_cast<int>(verts.size());
            verts.push_back(std::move(nvx));
        }
        cut.push_back(created);
        if (created != remap[in_old]) {
            Edge shortened;
            shortened.a = remap[in_old];
            shortened.b = created;
            shortened.p0 = e.p0;
            shortened.p1 = e.p1;
            if (!has_edge(edges, shortened.a, shortened.b))
                edges.push_back(shortened);
        }
    }

    for (int i = 0; i < nv; ++i)
        if (cls[i] == Side::On)
            mark_old(i);

    if (cut.size() < 3)
        return false;

    Vec3d center(0);
    for (int id : cut)
        center += verts[id].p;
    center /= static_cast<double>(cut.size());
    Vec3d u, v;
    plane_axes(pl.n, u, v);
    std::sort(cut.begin(), cut.end(), [&](int i, int j) {
        const Vec3d di = verts[i].p - center;
        const Vec3d dj = verts[j].p - center;
        const double ai = std::atan2(dot(di, v), dot(di, u));
        const double aj = std::atan2(dot(dj, v), dot(dj, u));
        if (ai != aj)
            return ai < aj;
        return sqr_length(di) < sqr_length(dj);
    });

    std::vector<int> ring;
    ring.reserve(cut.size());
    for (int id : cut) {
        if (!ring.empty() && (ring.back() == id || sqr_length(verts[id].p - verts[ring.back()].p) <= cell.tol * cell.tol))
            continue;
        ring.push_back(id);
    }
    if (ring.size() >= 2 && sqr_length(verts[ring.front()].p - verts[ring.back()].p) <= cell.tol * cell.tol)
        ring.pop_back();
    if (ring.size() < 3)
        return false;
    if (polygon_cross2(verts, ring) <= (cell.tol * cell.half) * (cell.tol * cell.half))
        return false;

    for (size_t k = 0; k < ring.size(); ++k) {
        const int i = ring[k];
        const int j = ring[(k + 1) % ring.size()];
        if (i == j || has_edge(edges, i, j))
            continue;
        Edge face;
        face.a = i;
        face.b = j;
        face.p0 = shared_plane(verts[i], verts[j], fresh, cell.planes, cell.tol);
        face.p1 = fresh;
        edges.push_back(face);
    }

    cell.planes.push_back(pl);
    cell.verts.swap(verts);
    cell.edges.swap(edges);
    return true;
}

bool positive_area(const Cell& cell, const std::vector<int>& ids)
{
    return polygon_cross2(cell.verts, ids) > (cell.tol * cell.half) * (cell.tol * cell.half);
}

void clip_neighborhood(Cell& cell, const std::vector<Vec3d>& points, int site, double radius,
                       const std::vector<std::pair<double, int>>& order)
{
    // Half-extent is the certificate radius, not the farthest neighbor.
    init_cube(cell, 0.5 * radius * (1.0 + 1e-8));
    for (const auto& item : order) {
        const Vec3d delta = points[item.second] - points[site];
        const double dist = std::sqrt(item.first);
        if (!(dist > 0))
            continue;
        Plane pl;
        pl.n = delta / dist;
        pl.d = -0.5 * dist;
        pl.src = item.second;
        clip_by_plane(cell, pl);
    }
}

std::vector<int> neighbors_from_cell(const Cell& cell)
{
    std::vector<int> nbrs;
    std::vector<int> on_plane;
    on_plane.reserve(8);
    for (int pi = 0; pi < static_cast<int>(cell.planes.size()); ++pi) {
        if (cell.planes[pi].src < 0)
            continue;
        on_plane.clear();
        for (int vi = 0; vi < static_cast<int>(cell.verts.size()); ++vi)
            if (std::abs(side_of(cell.planes[pi], cell.verts[vi].p)) <= cell.tol)
                on_plane.push_back(vi);
        if (positive_area(cell, on_plane))
            nbrs.push_back(cell.planes[pi].src);
    }
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    return nbrs;
}

std::vector<int> neighbors_of_site(const std::vector<Vec3d>& points, int site, double radius,
                                   const std::vector<std::pair<double, int>>& order)
{
    Cell cell;
    clip_neighborhood(cell, points, site, radius, order);
    return neighbors_from_cell(cell);
}

/// Cyclic order of real neighbors around plane `pi`, oriented so the face is
/// counterclockwise against the plane's outward normal. Empty when the face
/// walk does not close or does not meet every graph neighbor.
std::vector<int> real_neighbor_rotation(const Cell& cell, int pi, const std::vector<int>& graph_nbrs)
{
    const int nplanes = static_cast<int>(cell.planes.size());
    std::unordered_map<int, std::vector<std::pair<int, int>>> star;
    for (const Edge& e : cell.edges) {
        int other = -1;
        if (e.p0 == pi)
            other = e.p1;
        else if (e.p1 == pi)
            other = e.p0;
        else
            continue;
        if (other < 0 || other >= nplanes)
            continue;
        const int src = cell.planes[other].src;
        star[e.a].push_back({e.b, src});
        star[e.b].push_back({e.a, src});
    }
    if (star.size() < 3)
        return {};

    int start = std::numeric_limits<int>::max();
    for (const auto& entry : star)
        start = std::min(start, entry.first);

    std::vector<int> verts;
    std::vector<int> srcs;
    int prev = -1;
    int cur = start;
    bool returned = false;
    const int guard_max = static_cast<int>(star.size()) + 2;
    for (int step = 0; step < guard_max; ++step) {
        const auto it = star.find(cur);
        if (it == star.end())
            break;
        int next_v = -1;
        int next_src = -1;
        for (const auto& step_to : it->second) {
            if (step_to.first == prev)
                continue;
            if (next_v < 0 || step_to.first < next_v) {
                next_v = step_to.first;
                next_src = step_to.second;
            }
        }
        if (next_v < 0)
            break;
        verts.push_back(cur);
        srcs.push_back(next_src);
        prev = cur;
        cur = next_v;
        if (cur == start) {
            returned = true;
            break;
        }
    }
    if (!returned || verts.size() < 3)
        return {};

    Vec3d area(0);
    for (size_t i = 0; i < verts.size(); ++i) {
        const Vec3d& p = cell.verts[verts[i]].p;
        const Vec3d& q = cell.verts[verts[(i + 1) % verts.size()]].p;
        area += cross(p, q);
    }
    const bool outward_ccw = dot(area, cell.planes[pi].n) >= 0;

    std::vector<int> real;
    real.reserve(srcs.size());
    for (int src : srcs)
        if (src >= 0)
            real.push_back(src);
    if (!outward_ccw)
        std::reverse(real.begin(), real.end());

    std::vector<int> compact;
    compact.reserve(real.size());
    for (int src : real) {
        if (!std::binary_search(graph_nbrs.begin(), graph_nbrs.end(), src))
            continue;
        if (!compact.empty() && compact.back() == src)
            continue;
        compact.push_back(src);
    }
    if (compact.size() >= 2 && compact.front() == compact.back())
        compact.pop_back();

    std::vector<int> sorted_compact = compact;
    std::sort(sorted_compact.begin(), sorted_compact.end());
    sorted_compact.erase(std::unique(sorted_compact.begin(), sorted_compact.end()), sorted_compact.end());
    if (sorted_compact != graph_nbrs)
        return {};
    return compact;
}

DelaunayLink walk_path_or_cycle(const std::vector<int>& comp,
                                const std::unordered_map<int, std::vector<int>>& adj)
{
    int start = comp.front();
    bool has_end = false;
    for (int node : comp) {
        if (adj.at(node).size() <= 1 && (!has_end || node < start)) {
            start = node;
            has_end = true;
        }
    }

    DelaunayLink link;
    link.closed = !has_end && comp.size() >= 3;
    int prev = -1;
    int cur = start;
    link.neighbors.reserve(comp.size());
    for (size_t step = 0; step < comp.size(); ++step) {
        link.neighbors.push_back(cur);
        int next = -1;
        for (int nb : adj.at(cur)) {
            if (nb == prev)
                continue;
            next = nb;
            break;
        }
        if (next < 0) {
            link.closed = false;
            break;
        }
        if (next == start) {
            link.closed = link.neighbors.size() == comp.size() && adj.at(start).size() == 2;
            break;
        }
        prev = cur;
        cur = next;
    }
    if (link.neighbors.size() != comp.size())
        link.closed = false;
    return link;
}

std::vector<int> canonical_cycle(std::vector<int> face)
{
    if (face.empty())
        return face;
    const auto at = std::min_element(face.begin(), face.end());
    std::rotate(face.begin(), at, face.end());
    std::vector<int> reversed(face.size());
    reversed[0] = face[0];
    for (size_t i = 1; i < face.size(); ++i)
        reversed[i] = face[face.size() - i];
    if (reversed < face)
        return reversed;
    return face;
}

DelaunayLink outer_face(const std::vector<int>& comp,
                        const std::unordered_map<int, std::vector<int>>& adj,
                        const std::unordered_map<int, std::vector<int>>& rotation)
{
    auto next_after = [&](int at, int from) -> int {
        const std::vector<int>* rot = nullptr;
        const auto it = rotation.find(at);
        if (it != rotation.end() && !it->second.empty())
            rot = &it->second;
        else
            rot = &adj.at(at);
        for (size_t i = 0; i < rot->size(); ++i) {
            if ((*rot)[i] == from)
                return (*rot)[(i + 1) % rot->size()];
        }
        for (int nb : adj.at(at))
            if (nb != from)
                return nb;
        return -1;
    };

    std::set<std::pair<int, int>> unused;
    for (int u : comp) {
        for (int v : adj.at(u)) {
            if (u < v) {
                unused.emplace(u, v);
                unused.emplace(v, u);
            }
        }
    }

    DelaunayLink best;
    bool found = false;
    while (!unused.empty()) {
        const auto start_edge = *unused.begin();
        int from = start_edge.first;
        int to = start_edge.second;
        std::vector<int> face;
        bool closed_face = false;
        const int limit = static_cast<int>(unused.size()) + 2;
        for (int guard = 0; guard < limit; ++guard) {
            const std::pair<int, int> key{from, to};
            if (!unused.erase(key))
                break;
            face.push_back(from);
            const int next = next_after(to, from);
            if (next < 0 || next == to)
                break;
            from = to;
            to = next;
            if (from == start_edge.first && to == start_edge.second) {
                closed_face = true;
                break;
            }
        }
        if (!closed_face || face.size() < 2)
            continue;
        std::set<int> unique(face.begin(), face.end());
        if (unique.size() != face.size())
            continue;
        std::vector<int> canon = canonical_cycle(std::move(face));
        if (!found || canon.size() > best.neighbors.size()
            || (canon.size() == best.neighbors.size() && canon < best.neighbors)) {
            found = true;
            best.neighbors = std::move(canon);
            best.closed = best.neighbors.size() >= 3;
        }
    }
    return best;
}

std::vector<DelaunayLink> links_from_cell(const Cell& cell, const std::vector<int>& allowed)
{
    const int nplanes = static_cast<int>(cell.planes.size());
    std::unordered_map<int, std::vector<int>> adj;
    auto add_link = [&](int a, int b) {
        if (a < 0 || b < 0 || a == b)
            return;
        if (!std::binary_search(allowed.begin(), allowed.end(), a))
            return;
        if (!std::binary_search(allowed.begin(), allowed.end(), b))
            return;
        adj[a].push_back(b);
        adj[b].push_back(a);
    };
    for (const Edge& e : cell.edges) {
        if (e.p0 < 0 || e.p1 < 0 || e.p0 >= nplanes || e.p1 >= nplanes)
            continue;
        add_link(cell.planes[e.p0].src, cell.planes[e.p1].src);
    }
    for (auto& entry : adj) {
        auto& nbrs = entry.second;
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    if (adj.empty())
        return {};

    bool need_rotation = false;
    for (const auto& entry : adj) {
        if (entry.second.size() > 2)
            need_rotation = true;
    }
    std::unordered_map<int, int> src_plane;
    std::unordered_map<int, std::vector<int>> rotation;
    if (need_rotation) {
        for (int pi = 0; pi < nplanes; ++pi) {
            const int src = cell.planes[pi].src;
            if (src >= 0)
                src_plane.emplace(src, pi);
        }
        for (const auto& entry : adj) {
            const auto plane = src_plane.find(entry.first);
            if (plane == src_plane.end())
                continue;
            std::vector<int> rot = real_neighbor_rotation(cell, plane->second, entry.second);
            if (!rot.empty())
                rotation.emplace(entry.first, std::move(rot));
        }
    }

    std::vector<int> sources;
    sources.reserve(adj.size());
    for (const auto& entry : adj)
        sources.push_back(entry.first);
    std::sort(sources.begin(), sources.end());

    std::unordered_map<int, char> seen;
    std::vector<DelaunayLink> links;
    for (int src : sources) {
        if (seen[src])
            continue;
        std::vector<int> comp;
        std::vector<int> stack{src};
        seen[src] = 1;
        while (!stack.empty()) {
            const int node = stack.back();
            stack.pop_back();
            comp.push_back(node);
            for (int nb : adj[node]) {
                if (!seen[nb]) {
                    seen[nb] = 1;
                    stack.push_back(nb);
                }
            }
        }
        std::sort(comp.begin(), comp.end());

        int max_deg = 0;
        for (int node : comp)
            max_deg = std::max(max_deg, static_cast<int>(adj[node].size()));

        if (max_deg <= 2) {
            DelaunayLink link = walk_path_or_cycle(comp, adj);
            if (link.neighbors.size() >= 2)
                links.push_back(std::move(link));
            continue;
        }

        DelaunayLink link = outer_face(comp, adj, rotation);
        if (link.neighbors.size() >= 2) {
            links.push_back(std::move(link));
            continue;
        }
        for (int node : comp) {
            for (int nb : adj[node]) {
                if (node < nb)
                    links.push_back(DelaunayLink{{node, nb}, false});
            }
        }
    }
    return links;
}

DelaunayStar star_of_site(const std::vector<Vec3d>& points, int site, double radius,
                          const std::vector<std::pair<double, int>>& order)
{
    Cell cell;
    clip_neighborhood(cell, points, site, radius, order);
    DelaunayStar star;
    star.neighbors = neighbors_from_cell(cell);
    star.links = links_from_cell(cell, star.neighbors);
    return star;
}

template <class Fn>
void for_each_site_order(const std::vector<Vec3d>& points, double radius, Fn&& visit)
{
    const size_t n = points.size();
    if (n == 0 || !(radius > 0.0) || !std::isfinite(radius))
        return;
    if (n > static_cast<size_t>(std::numeric_limits<int>::max()))
        return;

    KDTree<Vec3d, int> tree;
    for (size_t i = 0; i < n; ++i)
        if (finite_point(points[i]))
            tree.insert(points[i], static_cast<int>(i));
    tree.build();

    const double query = std::nextafter(radius, std::numeric_limits<double>::infinity());
    const double limit2 = radius * radius * (1.0 + 1e-12);
    const double min2 = std::max(0.0, radius * radius * 1e-24);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0 || n < 64)
        threads = 1;
    threads = std::min(threads, static_cast<unsigned>(n));

    std::atomic<size_t> cursor{0};
    auto worker = [&]() {
        std::vector<std::pair<Vec3d, int>> hits;
        std::vector<std::pair<double, int>> order;
        while (true) {
            const size_t i = cursor.fetch_add(1, std::memory_order_relaxed);
            if (i >= n)
                break;
            if (!finite_point(points[i]))
                continue;
            hits.clear();
            tree.in_sphere(points[i], query, hits);
            order.clear();
            order.reserve(hits.size());
            for (const auto& hit : hits) {
                if (hit.second == static_cast<int>(i))
                    continue;
                const double d2 = sqr_length(points[i] - points[hit.second]);
                if (!(d2 > min2) || d2 > limit2 || !std::isfinite(d2))
                    continue;
                order.emplace_back(d2, hit.second);
            }
            if (order.empty())
                continue;
            std::sort(order.begin(), order.end());
            visit(static_cast<int>(i), order);
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned t = 0; t < threads; ++t)
        workers.emplace_back(worker);
    for (std::thread& worker_thread : workers)
        worker_thread.join();
}

/// Empty diametral ball: the midpoint is strictly inside every other bisector.
std::vector<int> gabriel_neighbors_of_site(const std::vector<Vec3d>& points, int site, double radius,
                                           const std::vector<std::pair<double, int>>& order)
{
    const double tol = std::max(1e-15, 1e-10 * 0.5 * radius);
    struct Bisector {
        Vec3d n;
        double d;
        int src;
    };
    std::vector<Bisector> planes;
    planes.reserve(order.size());
    for (const auto& item : order) {
        const Vec3d delta = points[item.second] - points[site];
        const double dist = std::sqrt(item.first);
        if (!(dist > 0) || !std::isfinite(dist))
            continue;
        planes.push_back({delta / dist, -0.5 * dist, item.second});
    }

    std::vector<int> nbrs;
    nbrs.reserve(planes.size());
    for (size_t i = 0; i < planes.size(); ++i) {
        const Vec3d midpoint = planes[i].n * (-planes[i].d);
        bool empty = true;
        for (size_t j = 0; j < planes.size(); ++j) {
            if (i == j)
                continue;
            if (dot(planes[j].n, midpoint) + planes[j].d >= -tol) {
                empty = false;
                break;
            }
        }
        if (empty)
            nbrs.push_back(planes[i].src);
    }
    return nbrs;
}

using SiteNeighbors = std::vector<int> (*)(const std::vector<Vec3d>&, int, double,
                                           const std::vector<std::pair<double, int>>&);

AMGraph3D radius_limited_edges(const std::vector<Vec3d>& points, double radius, SiteNeighbors select)
{
    AMGraph3D graph;
    const size_t n = points.size();
    for (size_t i = 0; i < n; ++i)
        graph.add_node(points[i]);
    if (select == nullptr)
        return graph;

    std::vector<std::vector<int>> adj(n);
    for_each_site_order(points, radius, [&](int site, const std::vector<std::pair<double, int>>& order) {
        adj[static_cast<size_t>(site)] = select(points, site, radius, order);
    });

    for (size_t i = 0; i < n; ++i)
        for (int j : adj[i])
            graph.connect_nodes(i, static_cast<AMGraph::NodeID>(j));
    return graph;
}

} // namespace

AMGraph3D delaunay_edges(const std::vector<CGLA::Vec3d>& points, double radius)
{
    return radius_limited_edges(points, radius, neighbors_of_site);
}

AMGraph3D gabriel_edges(const std::vector<CGLA::Vec3d>& points, double radius)
{
    return radius_limited_edges(points, radius, gabriel_neighbors_of_site);
}

std::vector<DelaunayStar> delaunay_stars(const std::vector<CGLA::Vec3d>& points, double radius)
{
    std::vector<DelaunayStar> stars(points.size());
    for_each_site_order(points, radius, [&](int site, const std::vector<std::pair<double, int>>& order) {
        stars[static_cast<size_t>(site)] = star_of_site(points, site, radius, order);
    });
    return stars;
}

} // namespace Geometry
