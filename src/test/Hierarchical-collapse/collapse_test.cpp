#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <GEL/Geometry/NeighborUtil.h>
#include <GEL/HMesh/HierarchicalReconstruction.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Util/ThreadPool.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
#include <vector>

using CGLA::Vec3d;
using HMesh::RSR::CollapseOpts;
using HMesh::RSR::Distance;
using HMesh::RSR::SimplifiedCloud;

namespace {

std::vector<Vec3d> unit_grid(const int side)
{
    std::vector<Vec3d> points;
    points.reserve(static_cast<size_t>(side) * side);
    for (int y = 0; y < side; ++y)
        for (int x = 0; x < side; ++x)
            points.emplace_back(x, y, 0.0);
    return points;
}

CollapseOpts options(const size_t iterations, const size_t reconstruction_neighbors)
{
    CollapseOpts opts;
    opts.max_iterations = iterations;
    opts.reduction_per_iteration = 0.5;
    opts.initial_neighbors = 5;
    opts.reconstruction_neighbors = reconstruction_neighbors;
    opts.distance = Distance::Euclidean;
    return opts;
}

struct GraphStats {
    size_t survivors = 0;
    size_t edges = 0;
    double max_length = 0;
};

GraphStats inspect(const SimplifiedCloud& cloud)
{
    GraphStats stats;
    stats.survivors = cloud.graph.positions.size();
    CHECK(cloud.cloud.points.size() == stats.survivors);
    CHECK(cloud.graph.normals.size() == stats.survivors);
    for (size_t i = 0; i < stats.survivors; ++i) {
        const Vec3d& point = cloud.graph.positions[i];
        CHECK(std::isfinite(point[0]));
        CHECK(std::isfinite(point[1]));
        CHECK(std::isfinite(point[2]));
        CHECK(cloud.graph.normals[i].length() == doctest::Approx(1.0).epsilon(1e-6));
        for (const size_t neighbor : cloud.graph.graph.neighbors(i)) {
            CHECK(neighbor != i);
            CHECK(neighbor < stats.survivors);
            if (i > neighbor)
                continue;
            ++stats.edges;
            const double length = (point - cloud.graph.positions[neighbor]).length();
            stats.max_length = std::max(stats.max_length, length);
        }
    }
    return stats;
}

SimplifiedCloud collapse_grid(const int side, const CollapseOpts& opts)
{
    const std::vector<Vec3d> points = unit_grid(side);
    const std::vector<Vec3d> normals(points.size(), Vec3d(0, 0, 1));
    return collapse_points(points, normals, opts).second;
}

Geometry::CloudNeighborhood search_grid(const std::vector<Vec3d>& points, const int neighbor_count)
{
    Geometry::Tree tree;
    std::vector<size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    Geometry::build_kd_tree_of_indices(points, indices, tree);
    Util::detail::ImmediatePool pool;
    return Geometry::search_neighborhoods(pool, points, tree, neighbor_count);
}

std::set<std::pair<size_t, size_t>> undirected_edges(const SimplifiedCloud& cloud)
{
    std::set<std::pair<size_t, size_t>> edges;
    for (size_t i = 0; i < cloud.graph.positions.size(); ++i) {
        for (const size_t neighbor : cloud.graph.graph.neighbors(i)) {
            if (i < neighbor)
                edges.emplace(i, neighbor);
        }
    }
    return edges;
}

void check_same_cloud(const SimplifiedCloud& first, const SimplifiedCloud& second)
{
    REQUIRE(first.graph.positions.size() == second.graph.positions.size());
    for (size_t i = 0; i < first.graph.positions.size(); ++i) {
        for (int axis = 0; axis < 3; ++axis)
            CHECK(first.graph.positions[i][axis] == doctest::Approx(second.graph.positions[i][axis]));
    }
    CHECK(undirected_edges(first) == undirected_edges(second));
}

}

TEST_CASE("ring-limited collapse keeps edges outside the one-ring")
{
    const CollapseOpts limited = options(2, 24);
    const CollapseOpts one_ring = options(2, 0);
    const GraphStats kept = inspect(collapse_grid(10, limited));
    const GraphStats local = inspect(collapse_grid(10, one_ring));

    INFO("kept survivors " << kept.survivors << " edges " << kept.edges << " max " << kept.max_length);
    INFO("one-ring survivors " << local.survivors << " edges " << local.edges << " max " << local.max_length);
    CHECK(kept.survivors < 100);
    CHECK(kept.survivors > 20);
    // A unit-grid one-ring is about length 1. Edges longer than that were kept
    // for reconstruction and must still be present after both iterations.
    CHECK(kept.max_length > 2.5);
    CHECK(kept.edges > local.edges * 2);
}

TEST_CASE("ring-limited collapse does not depend on hash order")
{
    const CollapseOpts opts = options(2, 24);
    const SimplifiedCloud first = collapse_grid(8, opts);
    const SimplifiedCloud second = collapse_grid(8, opts);
    REQUIRE(first.graph.positions.size() == second.graph.positions.size());
    for (size_t i = 0; i < first.graph.positions.size(); ++i) {
        for (int axis = 0; axis < 3; ++axis)
            CHECK(first.graph.positions[i][axis] == doctest::Approx(second.graph.positions[i][axis]));
    }
    CHECK(inspect(first).edges == inspect(second).edges);
}

TEST_CASE("one-ring collapse still removes vertices")
{
    const GraphStats stats = inspect(collapse_grid(8, options(1, 0)));
    CHECK(stats.survivors < 64);
    CHECK(stats.survivors > 16);
    CHECK(stats.edges > 0);
}

TEST_CASE("dense neighborhood collapse finishes")
{
    constexpr int side = 48;
    const std::vector<Vec3d> points = unit_grid(side);
    const std::vector<Vec3d> normals(points.size(), Vec3d(0, 0, 1));
    const auto elapsed_ms = [&](const CollapseOpts& opts) {
        const auto start = std::chrono::steady_clock::now();
        const SimplifiedCloud cloud = collapse_points(points, normals, opts).second;
        const auto stop = std::chrono::steady_clock::now();
        CHECK(cloud.graph.positions.size() < points.size());
        return std::chrono::duration<double, std::milli>(stop - start).count();
    };

    const double one_ring_ms = elapsed_ms(options(1, 0));
    const double limited_ms = elapsed_ms(options(1, 70));
    MESSAGE("collapse_points 48x48: one-ring " << one_ring_ms << " ms, k=70 " << limited_ms << " ms");
    CHECK(one_ring_ms > 0.0);
    CHECK(limited_ms > 0.0);
}

TEST_CASE("one-ring radius is the rank distance of a longer search")
{
    const std::vector<Vec3d> points = unit_grid(8);
    Geometry::CloudNeighborhood wide = search_grid(points, 24);
    Geometry::CloudNeighborhood ring = search_grid(points, 5);
    REQUIRE(wide.note_one_ring(5));
    CHECK_FALSE(ring.note_one_ring(24));
    REQUIRE(wide.points.size() == ring.points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        CHECK(wide.points[i].farthest_distance == doctest::Approx(wide.points[i].nearest.back().distance));
        CHECK(wide.points[i].one_ring_radius == doctest::Approx(ring.points[i].nearest.back().distance));
    }

    REQUIRE(wide.note_reconstruction_cap(24));
    const size_t cap_index = static_cast<size_t>(static_cast<double>(24) * (2.0 / 3.0));
    Geometry::CloudNeighborhood longer = search_grid(points, 30);
    REQUIRE(longer.note_reconstruction_cap(24));
    for (size_t i = 0; i < points.size(); ++i) {
        CHECK(wide.points[i].reconstruction_length_cap
              == doctest::Approx(wide.points[i].nearest[cap_index].distance));
        CHECK(longer.points[i].reconstruction_length_cap
              == doctest::Approx(wide.points[i].reconstruction_length_cap));
    }
}

TEST_CASE("collapse reuses an exact seed neighborhood")
{
    const std::vector<Vec3d> points = unit_grid(8);
    const std::vector<Vec3d> normals(points.size(), Vec3d(0, 0, 1));
    const CollapseOpts opts = options(2, 24);
    Geometry::CloudNeighborhood neighborhood = search_grid(points, 24);
    const SimplifiedCloud reused = collapse_points(points, normals, opts, neighborhood).second;
    const SimplifiedCloud fresh = collapse_points(points, normals, opts).second;
    check_same_cloud(reused, fresh);

    const Geometry::CloudNeighborhood ring = search_grid(points, 5);
    CHECK(neighborhood.one_ring_count == 5);
    for (size_t i = 0; i < points.size(); ++i)
        CHECK(neighborhood.points[i].one_ring_radius == doctest::Approx(ring.points[i].nearest.back().distance));
}

TEST_CASE("a longer neighborhood is not the collapse seed")
{
    const std::vector<Vec3d> points = unit_grid(8);
    const std::vector<Vec3d> normals(points.size(), Vec3d(0, 0, 1));
    const CollapseOpts opts = options(1, 0);
    Geometry::CloudNeighborhood neighborhood = search_grid(points, 24);
    const SimplifiedCloud ignored = collapse_points(points, normals, opts, neighborhood).second;
    const SimplifiedCloud fresh = collapse_points(points, normals, opts).second;
    check_same_cloud(ignored, fresh);
    CHECK(neighborhood.one_ring_count == 0);
    CHECK(neighborhood.searched_count == 24);
}

namespace HMesh::RSR {
size_t relax_reinserted_link(
    Manifold& manifold,
    VertexID kept,
    VertexID inserted,
    const CGLA::Vec3d& center,
    const CGLA::Vec3d& reference,
    std::vector<VertexID>* removed_vertices = nullptr);
}

namespace {

bool vertices_connected(const HMesh::Manifold& mesh, const HMesh::VertexID a, const HMesh::VertexID b)
{
    for (const HMesh::HalfEdgeID spoke : mesh.incident_halfedges(a)) {
        if (mesh.walker(spoke).vertex() == b)
            return true;
    }
    return false;
}

HMesh::VertexID vertex_for(const HMesh::Manifold& mesh, const HMesh::VertexAttributeVector<size_t>& origin, const size_t index)
{
    for (const HMesh::VertexID vertex : mesh.vertices()) {
        if (origin[vertex] == index)
            return vertex;
    }
    return HMesh::InvalidVertexID;
}

}

TEST_CASE("link repair flips an interior edge when the minimum angle rises")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0},
        {2, 0, 0},
        {1, -0.02, 0},
        {1, 2, 0},
    };
    const std::vector<size_t> faces = {0, 2, 1, 0, 1, 3};
    HMesh::Manifold mesh;
    const auto origin = HMesh::build_manifold(mesh, points, faces, 3);
    const HMesh::VertexID a = vertex_for(mesh, origin, 0);
    const HMesh::VertexID b = vertex_for(mesh, origin, 1);
    const HMesh::VertexID c = vertex_for(mesh, origin, 2);
    const HMesh::VertexID d = vertex_for(mesh, origin, 3);
    REQUIRE(vertices_connected(mesh, a, b));

    HMesh::RSR::relax_reinserted_link(mesh, c, d, Vec3d(1, 0, 0), Vec3d(0, 0, 1));

    CHECK_FALSE(vertices_connected(mesh, a, b));
    CHECK(vertices_connected(mesh, c, d));
    CHECK(mesh.no_faces() == 2);
}

TEST_CASE("link repair leaves an interior edge whose minimum angle does not rise")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
    };
    const std::vector<size_t> faces = {0, 1, 3, 0, 3, 2};
    HMesh::Manifold mesh;
    const auto origin = HMesh::build_manifold(mesh, points, faces, 3);
    const HMesh::VertexID a = vertex_for(mesh, origin, 0);
    const HMesh::VertexID b = vertex_for(mesh, origin, 1);
    const HMesh::VertexID c = vertex_for(mesh, origin, 2);
    const HMesh::VertexID d = vertex_for(mesh, origin, 3);
    REQUIRE(vertices_connected(mesh, a, d));

    HMesh::RSR::relax_reinserted_link(mesh, b, c, Vec3d(0.5, 0.5, 0), Vec3d(0, 0, 1));

    CHECK(vertices_connected(mesh, a, d));
    CHECK_FALSE(vertices_connected(mesh, b, c));
    CHECK(mesh.no_faces() == 2);
}

TEST_CASE("link repair leaves an inverted boundary edge in place")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0.5, 0.5, 0},
    };
    const std::vector<size_t> faces = {0, 1, 2, 0, 2, 3};
    HMesh::Manifold mesh;
    const auto origin = HMesh::build_manifold(mesh, points, faces, 3);
    const HMesh::VertexID center = vertex_for(mesh, origin, 0);
    const HMesh::VertexID right = vertex_for(mesh, origin, 2);
    const HMesh::VertexID inside = vertex_for(mesh, origin, 3);
    REQUIRE(mesh.no_faces() == 2);
    REQUIRE(mesh.valency(inside) == 2);

    std::vector<HMesh::VertexID> removed;
    HMesh::RSR::relax_reinserted_link(mesh, center, right, points[0], Vec3d(0, 0, 1), &removed);

    CHECK(mesh.no_faces() == 2);
    CHECK(mesh.no_vertices() == 4);
    CHECK(mesh.in_use(center));
    CHECK(mesh.in_use(right));
    CHECK(mesh.in_use(inside));
    CHECK(removed.empty());
}

TEST_CASE("link repair keeps a boundary edge that faces the reference normal")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {-0.5, 1, 0},
    };
    const std::vector<size_t> faces = {0, 1, 2, 0, 2, 3};
    HMesh::Manifold mesh;
    const auto origin = HMesh::build_manifold(mesh, points, faces, 3);
    const HMesh::VertexID center = vertex_for(mesh, origin, 0);
    const HMesh::VertexID right = vertex_for(mesh, origin, 2);

    HMesh::RSR::relax_reinserted_link(mesh, center, right, points[0], Vec3d(0, 0, 1));

    CHECK(mesh.no_faces() == 2);
}

TEST_CASE("link repair flips to a farther vertex when the minimum angle rises")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0},
        {-1, 0.01, 0},
        {1, 0.01, 0},
        {0, 5, 0},
        {0, -1, 0},
    };
    const std::vector<size_t> faces = {0, 2, 1, 1, 2, 3, 0, 1, 4};
    HMesh::Manifold mesh;
    const auto origin = HMesh::build_manifold(mesh, points, faces, 3);
    const HMesh::VertexID kept = vertex_for(mesh, origin, 0);
    const HMesh::VertexID left = vertex_for(mesh, origin, 1);
    const HMesh::VertexID right = vertex_for(mesh, origin, 2);
    const HMesh::VertexID far_vertex = vertex_for(mesh, origin, 3);
    const HMesh::VertexID inserted = vertex_for(mesh, origin, 4);
    REQUIRE(mesh.no_faces() == 3);
    REQUIRE(vertices_connected(mesh, left, right));
    REQUIRE_FALSE(vertices_connected(mesh, kept, far_vertex));

    const size_t flips = HMesh::RSR::relax_reinserted_link(mesh, kept, inserted, points[0], Vec3d(0, 0, 1));

    CHECK(vertices_connected(mesh, kept, far_vertex));
    CHECK_FALSE(vertices_connected(mesh, left, right));
    CHECK(mesh.no_faces() == 3);
    CHECK(mesh.no_vertices() == 5);
    CHECK(flips > 0);
    CHECK(flips <= 200);
}
