#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <GEL/Geometry/normal.h>
#include <GEL/HMesh/RSRExperimental.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <span>
#include <vector>

using CGLA::Vec3d;

namespace {

double angle_degrees(const Vec3d& normal, const Vec3d& reference)
{
    const double cosine = std::clamp(
        std::abs(CGLA::dot(CGLA::normalize(normal), CGLA::normalize(reference))), 0.0, 1.0);
    return std::acos(cosine) * 180.0 / std::numbers::pi;
}

const Vec3d z_axis(0, 0, 1);
const Vec3d x_axis(1, 0, 0);

}

TEST_CASE("clean plane keeps the smallest neighborhood")
{
    std::vector<Vec3d> points = {
        {0, 0, 0},   {1.0, 0.1, 0}, {0.2, 1.4, 0}, {-0.7, 0.4, 0},
        {1.6, -0.3, 0}, {-1.1, -0.6, 0}, {0.4, -1.2, 0}, {2.2, 0.8, 0},
    };
    for (int i = 0; i < 24; ++i)
        points.emplace_back((i % 5) * 0.37 + 2.5, (i / 5) * 0.53 - 1.4, 0.0);

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(points);

    CHECK(estimated.k == Geometry::ADAPTIVE_NORMAL_SCHEDULE.front());
    CHECK(angle_degrees(estimated.normal, z_axis) < 1.0);
}

TEST_CASE("noise prefers a larger neighborhood than the smallest prefix")
{
    const double c = 0.5;
    const double s = std::sqrt(3.0) / 2.0;
    auto tilt = [c, s](double x, double y, double z) {
        return Vec3d(x, y * c - z * s, y * s + z * c);
    };

    std::vector<Vec3d> points;
    for (const double x : {-0.30, 0.30})
        for (const double y : {-0.25, 0.25})
            for (const double z : {-0.20, 0.20})
                points.push_back(tilt(x, y, z));
    REQUIRE(points.size() == 8);

    for (int ix = -5; ix <= 5; ++ix)
        for (int iy = -5; iy <= 5; ++iy)
            points.emplace_back(ix * 0.4, iy * 0.4, 0.0);

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(points);
    const Geometry::AdaptiveNormal small = Geometry::estimateNormalAdaptive(std::span(points.data(), 8));

    CHECK(estimated.k > Geometry::ADAPTIVE_NORMAL_SCHEDULE.front());
    CHECK(angle_degrees(small.normal, z_axis) > 40.0);
    CHECK(angle_degrees(estimated.normal, z_axis) < 15.0);
    CHECK(angle_degrees(estimated.normal, z_axis) < angle_degrees(small.normal, z_axis));
}

TEST_CASE("clean curved patch stays on the smallest neighborhood")
{
    std::vector<Vec3d> points;
    for (int ix = -8; ix <= 8; ++ix) {
        for (int iy = -8; iy <= 8; ++iy) {
            const double x = ix * 0.15;
            const double y = iy * 0.15;
            points.emplace_back(x, y, 1.5 * (x * x + y * y));
        }
    }
    std::ranges::sort(points, [](const Vec3d& a, const Vec3d& b) {
        return CGLA::sqr_length(a) < CGLA::sqr_length(b);
    });

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(points);

    CHECK(estimated.k == Geometry::ADAPTIVE_NORMAL_SCHEDULE.front());
    CHECK(angle_degrees(estimated.normal, z_axis) < 5.0);
}

TEST_CASE("a flat strip grows toward a rounder neighborhood")
{
    std::vector<Vec3d> points = {
        {0, 0, 0},
        {0.20, 0.10, 0},
        {-0.20, -0.10, 0},
        {0.36, -0.14, 0},
        {-0.36, 0.14, 0},
        {0.12, 0.18, 0},
        {-0.12, -0.18, 0},
        {0.40, 0.04, 0},
    };
    for (int ix = -5; ix <= 5; ++ix) {
        for (int iy = -5; iy <= 5; ++iy) {
            const Vec3d point(ix * 0.32, iy * 0.32, 0.0);
            if (point.length() >= 0.50)
                points.push_back(point);
        }
    }
    std::ranges::sort(points, [](const Vec3d& a, const Vec3d& b) {
        return CGLA::sqr_length(a) < CGLA::sqr_length(b);
    });

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(points);

    CHECK(estimated.k >= 20);
    CHECK(angle_degrees(estimated.normal, z_axis) < 5.0);
}

TEST_CASE("a line is rejected and the crease keeps the sheet normal")
{
    std::vector<Vec3d> line;
    for (int i = 0; i < 30; ++i)
        line.emplace_back((i - 15) * 0.07, 0.0, 0.0);

    const Geometry::AdaptiveNormal along_line = Geometry::estimateNormalAdaptive(line);
    CHECK(along_line.k == Geometry::ADAPTIVE_NORMAL_SCHEDULE.front());
    CHECK(angle_degrees(along_line.normal, x_axis) > 80.0);

    std::vector<Vec3d> crease;
    for (int i = 0; i < 27; ++i)
        crease.emplace_back((i - 13) * 0.08, 0.0, 0.0);
    for (int i = 0; i < 13; ++i)
        crease.emplace_back((i - 6) * 0.15, (i % 2 == 0) ? 1.0 : -1.0, 0.0);

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(crease);
    CHECK(estimated.k >= 30);
    CHECK(angle_degrees(estimated.normal, z_axis) < 5.0);
}

TEST_CASE("orientation weight penalizes an edge along the normal")
{
    const Vec3d up(0, 0, 1);
    const Vec3d down(0, 0, -1);
    const Vec3d side(1, 0, 0);
    const double scale = 2.0;

    const double tangent = Geometry::orientation_edge_weight(up, up, Vec3d(scale, 0, 0), scale);
    const double tunnel = Geometry::orientation_edge_weight(up, down, Vec3d(0, 0, scale), scale);
    const double corner = Geometry::orientation_edge_weight(up, side, Vec3d(0, scale, 0), scale);
    const double unscaled = Geometry::orientation_edge_weight(up, side, Vec3d(0, 0, scale), 0.0);

    CHECK(tangent < 1e-8);
    CHECK(std::abs(tunnel - 2.0) < 1e-8);
    CHECK(tunnel > corner);
    // No neighborhood radius: the tunnel term is left out and only the angle remains.
    CHECK(std::abs(unscaled - 1.0) < 1e-8);
}

TEST_CASE("orientation cone contains the sheet in front of a tunnel")
{
    const Vec3d up(0, 0, 1);
    const Vec3d parent_offset(0, 0, 0.15);
    CHECK(Geometry::orientation_edge_is_tough(up, up, parent_offset));
    CHECK_FALSE(Geometry::orientation_edge_is_tough(up, up, Vec3d(1, 0, 0)));

    const Vec3d axis = Geometry::orientation_cone_axis(-up, parent_offset);
    CHECK(CGLA::dot(axis, up) > 0.9);
    CHECK(Geometry::orientation_cone_contains(parent_offset, axis, 0.15));
    CHECK(Geometry::orientation_cone_contains(Vec3d(0.04, 0.0, 0.15), axis, 0.15));
    CHECK_FALSE(Geometry::orientation_cone_contains(Vec3d(0.05, 0.0, 0.0), axis, 0.15));
    CHECK_FALSE(Geometry::orientation_cone_contains(Vec3d(0.0, 0.0, 0.5), axis, 0.15));
}

TEST_CASE("orientation cluster score favors tangent and tunnel evidence")
{
    const Vec3d up(0, 0, 1);
    const Vec3d down(0, 0, -1);
    const Vec3d tangent_edge(2, 0, 0);
    const Vec3d tunnel_edge(0, 0, 3);

    const auto tangent_parallel = Geometry::orientation_cluster_evidence(up, up, tangent_edge);
    const auto tangent_anti = Geometry::orientation_cluster_evidence(up, down, tangent_edge);
    const auto tunnel_parallel = Geometry::orientation_cluster_evidence(up, up, tunnel_edge);
    const auto tunnel_anti = Geometry::orientation_cluster_evidence(up, down, tunnel_edge);
    const auto angled = Geometry::orientation_cluster_evidence(up, up, Vec3d(1, 0, 1));
    const auto degenerate = Geometry::orientation_cluster_evidence(up, up, Vec3d(0, 0, 0));
    const auto zero_normal = Geometry::orientation_cluster_evidence(Vec3d(0, 0, 0), up, tangent_edge);

    CHECK(std::abs(tangent_parallel.score - 1.0) < 1e-8);
    CHECK(std::abs(tangent_anti.score - 1.0) < 1e-8);
    CHECK(std::abs(tunnel_parallel.score - 1.0) < 1e-8);
    CHECK(std::abs(tunnel_anti.score - 1.0) < 1e-8);
    CHECK(std::abs(angled.score - 0.5) < 1e-8);
    CHECK(degenerate.score == 0.0);
    CHECK(degenerate.signed_evidence == 0.0);
    CHECK(zero_normal.score == 0.0);
    CHECK(zero_normal.signed_evidence == 0.0);

    CHECK(std::abs(tangent_parallel.signed_evidence - 1.0) < 1e-8);
    CHECK(std::abs(tangent_anti.signed_evidence - (-1.0)) < 1e-8);
    CHECK(std::abs(tunnel_parallel.signed_evidence - (-1.0)) < 1e-8);
    CHECK(std::abs(tunnel_anti.signed_evidence - 1.0) < 1e-8);
}

TEST_CASE("two close sheets stay coherent and agree")
{
    std::vector<Vec3d> points;
    points.reserve(72);
    for (int sheet = 0; sheet < 2; ++sheet) {
        const double z = sheet * 0.15;
        for (int y = 0; y < 6; ++y)
            for (int x = 0; x < 6; ++x)
                points.emplace_back(x * 0.5, y * 0.5, z);
    }

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 20;
    opts.cone_check = false;

    std::vector<Vec3d> normals;
    HMesh::RSR::point_cloud_normals(points, opts, normals);
    REQUIRE(normals.size() == points.size());

    const auto mean_of = [&](int begin, int end) {
        Vec3d sum(0);
        for (int i = begin; i < end; ++i)
            sum += normals[static_cast<size_t>(i)];
        return CGLA::normalize(sum);
    };
    const Vec3d mean0 = mean_of(0, 36);
    const Vec3d mean1 = mean_of(36, 72);

    for (int i = 0; i < 36; ++i)
        CHECK(CGLA::dot(CGLA::normalize(normals[static_cast<size_t>(i)]), mean0) > 0.9);
    for (int i = 36; i < 72; ++i)
        CHECK(CGLA::dot(CGLA::normalize(normals[static_cast<size_t>(i)]), mean1) > 0.9);
    // The gap is a tunnel. Joining it must not reverse one sheet.
    CHECK(CGLA::dot(mean0, mean1) > 0.5);
}

TEST_CASE("a wing over a body keeps the body's orientation")
{
    const double spacing = 0.08;
    std::vector<Vec3d> points;
    const auto add_grid = [&](int x0, int x1, int y0, int y1, double z) {
        for (int ix = x0; ix <= x1; ++ix)
            for (int iy = y0; iy <= y1; ++iy)
                points.emplace_back(ix * spacing, iy * spacing, z);
    };
    add_grid(0, 25, 0, 20, 0.0);
    const int body_count = static_cast<int>(points.size());
    add_grid(5, 25, 4, 16, 0.12);
    const int wing_end = static_cast<int>(points.size());
    points.emplace_back(0.4, 0.7, 0.04);
    points.emplace_back(0.4, 0.7, 0.08);

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 25;
    opts.cone_check = false;

    std::vector<Vec3d> normals;
    HMesh::RSR::point_cloud_normals(points, opts, normals);
    REQUIRE(normals.size() == points.size());

    const auto mean_of = [&](int begin, int end) {
        Vec3d sum(0);
        for (int i = begin; i < end; ++i)
            sum += normals[static_cast<size_t>(i)];
        return CGLA::normalize(sum);
    };
    const Vec3d body = mean_of(0, body_count);
    const Vec3d wing = mean_of(body_count, wing_end);
    for (int i = 0; i < body_count; ++i)
        CHECK(CGLA::dot(CGLA::normalize(normals[static_cast<size_t>(i)]), body) > 0.8);
    for (int i = body_count; i < wing_end; ++i)
        CHECK(CGLA::dot(CGLA::normalize(normals[static_cast<size_t>(i)]), wing) > 0.8);
    CHECK(CGLA::dot(body, wing) > 0.5);
}

TEST_CASE("cone check faces a dense second sheet away from its parent")
{
    std::vector<Vec3d> points;
    points.reserve(128);
    for (int sheet = 0; sheet < 2; ++sheet) {
        const double z = sheet * 0.15;
        for (int y = 0; y < 8; ++y)
            for (int x = 0; x < 8; ++x)
                points.emplace_back(x * 0.05, y * 0.05, z);
    }

    const auto mean_of = [](const std::vector<Vec3d>& normals, int begin, int end) {
        Vec3d sum(0);
        for (int i = begin; i < end; ++i)
            sum += normals[static_cast<size_t>(i)];
        return CGLA::normalize(sum);
    };
    const auto coherent = [](const std::vector<Vec3d>& normals, int begin, int end, const Vec3d& mean) {
        for (int i = begin; i < end; ++i)
            CHECK(CGLA::dot(CGLA::normalize(normals[static_cast<size_t>(i)]), mean) > 0.9);
    };

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 48;
    opts.cone_check = true;

    std::vector<Vec3d> faced;
    std::vector<HMesh::RSR::OrientationTreeEdge> tree;
    HMesh::RSR::point_cloud_normals(
        points, opts, faced, HMesh::RSR::NormalOrientation::MinimumSpanningTree, &tree);
    REQUIRE(faced.size() == points.size());
    CHECK(tree.size() + 1 == points.size());
    const int invoked = static_cast<int>(std::ranges::count_if(
        tree, [](const HMesh::RSR::OrientationTreeEdge& edge) { return edge.cone_invoked; }));
    CHECK(invoked >= 1);
    const Vec3d faced0 = mean_of(faced, 0, 64);
    const Vec3d faced1 = mean_of(faced, 64, 128);
    coherent(faced, 0, 64, faced0);
    coherent(faced, 64, 128, faced1);
    CHECK(CGLA::dot(faced0, faced1) < -0.5);

    opts.cone_check = false;
    std::vector<Vec3d> agreed;
    HMesh::RSR::point_cloud_normals(points, opts, agreed);
    REQUIRE(agreed.size() == points.size());
    const Vec3d agreed0 = mean_of(agreed, 0, 64);
    const Vec3d agreed1 = mean_of(agreed, 64, 128);
    coherent(agreed, 0, 64, agreed0);
    coherent(agreed, 64, 128, agreed1);
    CHECK(CGLA::dot(agreed0, agreed1) > 0.5);
}

TEST_CASE("mst orientation agrees across a plane")
{
    std::vector<Vec3d> points;
    points.reserve(64);
    for (int y = 0; y < 8; ++y)
        for (int x = 0; x < 8; ++x)
            points.emplace_back(x * 0.25, y * 0.25, 0.0);

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 12;

    std::vector<Vec3d> normals;
    std::vector<HMesh::RSR::OrientationTreeEdge> tree;
    HMesh::RSR::point_cloud_normals(
        points, opts, normals, HMesh::RSR::NormalOrientation::MinimumSpanningTree, &tree);
    REQUIRE(normals.size() == points.size());
    CHECK(tree.size() + 1 == points.size());
    CHECK(std::ranges::none_of(
        tree, [](const HMesh::RSR::OrientationTreeEdge& edge) { return edge.cone_invoked; }));

    Vec3d sum(0);
    for (const auto& normal : normals) {
        CHECK(std::abs(normal.length() - 1.0) < 1e-6);
        sum += normal;
    }
    const Vec3d mean = CGLA::normalize(sum);
    for (const auto& normal : normals)
        CHECK(CGLA::dot(normal, mean) > 0.9);
}

TEST_CASE("spt orientation agrees across a plane")
{
    std::vector<Vec3d> points;
    points.reserve(64);
    for (int y = 0; y < 8; ++y)
        for (int x = 0; x < 8; ++x)
            points.emplace_back(x * 0.25, y * 0.25, 0.0);

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 12;

    std::vector<Vec3d> mst_normals;
    std::vector<Vec3d> spt_normals;
    std::vector<HMesh::RSR::OrientationTreeEdge> tree;
    HMesh::RSR::point_cloud_normals(
        points, opts, mst_normals, HMesh::RSR::NormalOrientation::MinimumSpanningTree);
    HMesh::RSR::point_cloud_normals(
        points, opts, spt_normals, HMesh::RSR::NormalOrientation::ShortestPathTree, &tree);
    REQUIRE(spt_normals.size() == points.size());
    CHECK(tree.size() + 1 == points.size());
    CHECK(std::ranges::none_of(
        tree, [](const HMesh::RSR::OrientationTreeEdge& edge) { return edge.cone_invoked; }));

    Vec3d sum(0);
    for (size_t i = 0; i < spt_normals.size(); ++i) {
        CHECK(std::abs(spt_normals[i].length() - 1.0) < 1e-6);
        CHECK(CGLA::dot(spt_normals[i], mst_normals[i]) > 0.999);
        sum += spt_normals[i];
    }
    const Vec3d mean = CGLA::normalize(sum);
    for (const auto& normal : spt_normals)
        CHECK(CGLA::dot(normal, mean) > 0.9);
}

TEST_CASE("mst joins two sheets with one tunnel edge")
{
    std::vector<Vec3d> points;
    points.reserve(72);
    for (int sheet = 0; sheet < 2; ++sheet) {
        const double z = sheet * 0.15;
        for (int y = 0; y < 6; ++y)
            for (int x = 0; x < 6; ++x)
                points.emplace_back(x * 0.5, y * 0.5, z);
    }

    HMesh::RSR::RSROpts opts;
    opts.dist = HMesh::RSR::Distance::Euclidean;
    opts.num_neighbors = 20;
    opts.cone_check = false;

    std::vector<Vec3d> normals;
    std::vector<HMesh::RSR::OrientationTreeEdge> tree;
    HMesh::RSR::point_cloud_normals(
        points, opts, normals, HMesh::RSR::NormalOrientation::MinimumSpanningTree, &tree);
    REQUIRE(tree.size() + 1 == points.size());

    int crossings = 0;
    for (const auto& edge : tree) {
        const bool parent_top = edge.parent >= 36;
        const bool child_top = edge.child >= 36;
        if (parent_top != child_top)
            ++crossings;
    }
    // Tangent edges cost nothing. One tunnel joins the sheets. The rest stay out.
    CHECK(crossings == 1);
}

TEST_CASE("fewer than eight points uses the whole cloud")
{
    const std::vector<Vec3d> points = {
        {0, 0, 0}, {1, 0, 0}, {0, 2, 0}, {0.3, 0.4, 0}, {-0.5, 0.7, 0},
    };

    const Geometry::AdaptiveNormal estimated = Geometry::estimateNormalAdaptive(points);

    CHECK(estimated.k == 5);
    CHECK(angle_degrees(estimated.normal, z_axis) < 1.0);
}
