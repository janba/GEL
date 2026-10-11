#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <GEL/Geometry/normal.h>

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
    CHECK(estimated.k == 40);
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
    CHECK(std::abs(unscaled - 1.0) < 1e-8);
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
