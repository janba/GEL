#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/curvature.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

using CGLA::Mat2x2d;
using CGLA::Mat3x3d;
using CGLA::Mat4x4d;
using CGLA::MatN;
using CGLA::Vec2d;
using CGLA::Vec3d;
using CGLA::Vec3f;
using CGLA::Vec4d;
using CGLA::VecN;
using CGLA::cross;
using CGLA::dot;
using CGLA::length;
using CGLA::normalize;

namespace {

double alignment(const Vec3d& a, const Vec3d& b)
{
    return std::abs(dot(normalize(a), normalize(b)));
}

template <class MT>
std::vector<double> stored_eigenvalues(const MT& L, int count)
{
    std::vector<double> values;
    values.reserve(count);
    for (int i = 0; i < count; ++i)
        values.push_back(L[i][i]);
    std::sort(values.begin(), values.end());
    return values;
}

template <class MT>
void check_eigenpairs(const MT& A, const MT& Q, const MT& L, int count, double tolerance)
{
    for (int i = 0; i < count; ++i) {
        const auto residual = A * Q[i] - L[i][i] * Q[i];
        CHECK(length(residual) < tolerance * (1.0 + std::abs(L[i][i])));
        CHECK(std::abs(length(Q[i]) - 1.0) < tolerance);
        for (int j = i + 1; j < count; ++j)
            CHECK(std::abs(dot(Q[i], Q[j])) < tolerance);
    }
}

void closed_form_2x2(const Mat2x2d& A, double& smaller, double& larger)
{
    const double a = A[0][0];
    const double off = 0.5 * (A[0][1] + A[1][0]);
    const double c = A[1][1];
    const double disc = std::sqrt((a - c) * (a - c) + 4.0 * off * off);
    larger = 0.5 * (a + c + disc);
    smaller = 0.5 * (a + c - disc);
}

Vec3d paraboloid_point(double x, double y, double cu, double cuv, double cv,
                        const Vec3d& origin, const Vec3d& axis_x, const Vec3d& axis_y, const Vec3d& axis_n)
{
    const double height = 0.5 * cu * x * x + cuv * x * y + 0.5 * cv * y * y;
    return origin + x * axis_x + y * axis_y + height * axis_n;
}

/// outer_count link edges, starting at the first, get a triangle on the far side.
/// The returned ids are those outer vertices, in edge order. A short count leaves
/// the remaining opposite vertices missing.
HMesh::VertexID make_paraboloid(HMesh::Manifold& mesh,
                                const std::vector<double>& angles,
                                const std::vector<double>& radii,
                                double cu, double cuv, double cv,
                                const Vec3d& origin,
                                const Vec3d& axis_x,
                                const Vec3d& axis_y,
                                const Vec3d& axis_n,
                                int outer_count = 0,
                                std::vector<HMesh::VertexID>* outers = nullptr)
{
    const int count = static_cast<int>(angles.size());
    std::vector<Vec3d> vertices;
    vertices.reserve(count + 1 + std::max(outer_count, 0));
    vertices.push_back(origin);
    for (int i = 0; i < count; ++i) {
        const double x = radii[i] * std::cos(angles[i]);
        const double y = radii[i] * std::sin(angles[i]);
        vertices.push_back(paraboloid_point(x, y, cu, cuv, cv, origin, axis_x, axis_y, axis_n));
    }

    std::vector<std::size_t> faces;
    faces.reserve(3 * (count + std::max(outer_count, 0)));
    for (int i = 0; i < count; ++i) {
        faces.push_back(0);
        faces.push_back(1 + i);
        faces.push_back(1 + (i + 1) % count);
    }
    const int built_outers = std::max(0, std::min(outer_count, count));
    for (int i = 0; i < built_outers; ++i) {
        const int i0 = 1 + i;
        const int i1 = 1 + (i + 1) % count;
        const Vec3d local0 = vertices[i0] - origin;
        const Vec3d local1 = vertices[i1] - origin;
        const double x0 = dot(local0, axis_x);
        const double y0 = dot(local0, axis_y);
        const double x1 = dot(local1, axis_x);
        const double y1 = dot(local1, axis_y);
        const double mx = 0.5 * (x0 + x1);
        const double my = 0.5 * (y0 + y1);
        const double mid = std::sqrt(mx * mx + my * my);
        const double scale = mid > 1e-8 ? 1.0 + 0.12 / mid : 1.0;
        vertices.push_back(paraboloid_point(scale * mx, scale * my, cu, cuv, cv, origin, axis_x, axis_y, axis_n));
        const std::size_t outer = vertices.size() - 1;
        // Opposite winding to the inner face so the shared edge is traversed both ways
        // and the normal stays aligned with axis_n.
        faces.push_back(i1);
        faces.push_back(i0);
        faces.push_back(outer);
    }

    const auto original_index = HMesh::build_manifold(mesh, vertices, faces, 3);
    HMesh::VertexID center = HMesh::InvalidVertexID;
    if (outers != nullptr)
        outers->clear();
    for (const auto vertex : mesh.vertices()) {
        const std::size_t original = original_index[vertex];
        if (original == 0)
            center = vertex;
        else if (outers != nullptr && original >= static_cast<std::size_t>(count + 1))
            outers->push_back(vertex);
    }
    return center;
}

void check_principal(const HMesh::PrincipalCurvatures& pc,
                     double expected_min, double expected_max,
                     const Vec3d& expected_min_dir, const Vec3d& expected_max_dir,
                     double curvature_tol, double direction_cos)
{
    CHECK(std::abs(pc.min_curvature - expected_min) < curvature_tol);
    CHECK(std::abs(pc.max_curvature - expected_max) < curvature_tol);
    CHECK(alignment(pc.min_curv_direction, expected_min_dir) > direction_cos);
    CHECK(alignment(pc.max_curv_direction, expected_max_dir) > direction_cos);
    CHECK(std::abs(dot(pc.min_curv_direction, pc.max_curv_direction)) < 1e-6);
    CHECK(std::abs(length(pc.min_curv_direction) - 1.0) < 1e-6);
    CHECK(std::abs(length(pc.max_curv_direction) - 1.0) < 1e-6);
}

}

TEST_CASE("lin_solve recovers a five-dimensional solution")
{
    MatN<5> A(0);
    VecN<5> expected(0);
    for (unsigned int i = 0; i < 5; ++i) {
        expected[i] = 0.5 * static_cast<double>(i) - 1.25;
        for (unsigned int j = 0; j < 5; ++j)
            A[i][j] = (i == j) ? 3.0 + 0.2 * i : 0.05 * (static_cast<int>(i) - static_cast<int>(j));
    }
    // A zero on the first pivot forces a row exchange.
    A[0][0] = 0;

    VecN<5> b(0);
    for (unsigned int i = 0; i < 5; ++i)
        for (unsigned int j = 0; j < 5; ++j)
            b[i] += A[i][j] * expected[j];

    const VecN<5> solved = CGLA::lin_solve(A, b);
    for (unsigned int i = 0; i < 5; ++i)
        CHECK(solved[i] == doctest::Approx(expected[i]).epsilon(1e-9));
}

TEST_CASE("ls_solve recovers three and five coefficients")
{
    const std::vector<Vec3d> rows3 = {
        Vec3d(1, 0, 0), Vec3d(0, 1, 0), Vec3d(0, 0, 1), Vec3d(1, 1, 1)};
    const std::vector<double> rhs3 = {2.0, -1.0, 4.0, 5.0};
    const Vec3d solved3 = CGLA::ls_solve(rows3, rhs3);
    CHECK(solved3[0] == doctest::Approx(2.0).epsilon(1e-8));
    CHECK(solved3[1] == doctest::Approx(-1.0).epsilon(1e-8));
    CHECK(solved3[2] == doctest::Approx(4.0).epsilon(1e-8));

    VecN<5> expected(0);
    expected[0] = 0.2;
    expected[1] = -0.4;
    expected[2] = 1.5;
    expected[3] = -0.3;
    expected[4] = 0.8;

    std::vector<VecN<5>> rows;
    std::vector<double> rhs;
    const double samples[8][2] = {
        {0.3, 0.1}, {-0.2, 0.4}, {0.5, -0.3}, {-0.4, -0.2},
        {0.1, 0.6}, {-0.6, 0.2}, {0.25, -0.55}, {-0.15, -0.45}};
    for (const auto& sample : samples) {
        const double x = sample[0];
        const double y = sample[1];
        VecN<5> row(0);
        row[0] = x;
        row[1] = y;
        row[2] = 0.5 * x * x;
        row[3] = x * y;
        row[4] = 0.5 * y * y;
        double height = 0;
        for (unsigned int i = 0; i < 5; ++i)
            height += row[i] * expected[i];
        rows.push_back(row);
        rhs.push_back(height);
    }

    const VecN<5> solved = CGLA::ls_solve(rows, rhs);
    for (unsigned int i = 0; i < 5; ++i)
        CHECK(solved[i] == doctest::Approx(expected[i]).epsilon(1e-8));
}

TEST_CASE("lin_solve and ls_solve accept the built-in vectors")
{
    const Mat2x2d rows2(0, 2, 3, 0);
    const Vec2d solved2 = CGLA::lin_solve(rows2, Vec2d(4, 9));
    CHECK(solved2[0] == doctest::Approx(3.0).epsilon(1e-9));
    CHECK(solved2[1] == doctest::Approx(2.0).epsilon(1e-9));

    Mat3x3d rows3(0);
    rows3[0] = Vec3d(0, 1, 0);
    rows3[1] = Vec3d(0, 0, 2);
    rows3[2] = Vec3d(4, 0, 0);
    const Vec3d solved3 = CGLA::lin_solve(rows3, Vec3d(3, 8, 12));
    CHECK(solved3[0] == doctest::Approx(3.0).epsilon(1e-9));
    CHECK(solved3[1] == doctest::Approx(3.0).epsilon(1e-9));
    CHECK(solved3[2] == doctest::Approx(4.0).epsilon(1e-9));

    Mat4x4d rows4(0);
    const Vec4d expected4(1.0, -2.0, 0.5, 4.0);
    for (unsigned int i = 0; i < 4; ++i) {
        rows4[i][i] = 2.0 + 0.1 * i;
        rows4[i][(i + 1) % 4] = -0.25;
    }
    rows4[0][0] = 0;
    Vec4d rhs4(0);
    for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 4; ++j)
            rhs4[i] += rows4[i][j] * expected4[j];
    const Vec4d solved4 = CGLA::lin_solve(rows4, rhs4);
    for (unsigned int i = 0; i < 4; ++i)
        CHECK(solved4[i] == doctest::Approx(expected4[i]).epsilon(1e-9));

    const std::vector<Vec3f> samples = {
        Vec3f(1, 0, 0), Vec3f(0, 1, 0), Vec3f(0, 0, 1), Vec3f(1, 1, 1)};
    const std::vector<float> rhs = {2.f, -1.f, 4.f, 5.f};
    const Vec3f solved_f = CGLA::ls_solve(samples, rhs);
    CHECK(solved_f[0] == doctest::Approx(2.f).epsilon(1e-5));
    CHECK(solved_f[1] == doctest::Approx(-1.f).epsilon(1e-5));
    CHECK(solved_f[2] == doctest::Approx(4.f).epsilon(1e-5));
}

TEST_CASE("power_eigensolution separates equal opposite eigenvalues")
{
    const double cosine = std::cos(0.4);
    const double sine = std::sin(0.4);
    const Mat2x2d turn(cosine, -sine, sine, cosine);
    const Mat2x2d diagonal(2.0, 0.0, 0.0, -2.0);
    const Mat2x2d A = turn * diagonal * CGLA::transpose(turn);

    Mat2x2d Q, L;
    const int count = CGLA::power_eigensolution(A, Q, L);
    REQUIRE(count == 2);
    check_eigenpairs(A, Q, L, count, 1e-5);

    const auto values = stored_eigenvalues(L, count);
    CHECK(values[0] == doctest::Approx(-2.0).epsilon(1e-5));
    CHECK(values[1] == doctest::Approx(2.0).epsilon(1e-5));
}

TEST_CASE("power_eigensolution reports nonzero rank")
{
    Mat3x3d positive(0);
    positive[0][0] = 3;
    positive[1][1] = 1;
    Mat3x3d Q, L;
    int count = CGLA::power_eigensolution(positive, Q, L);
    REQUIRE(count == 2);
    check_eigenpairs(positive, Q, L, count, 1e-5);
    const auto positive_values = stored_eigenvalues(L, count);
    CHECK(positive_values[0] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(positive_values[1] == doctest::Approx(3.0).epsilon(1e-5));

    Mat3x3d negative(0);
    negative[0][0] = -2;
    negative[1][1] = -1;
    count = CGLA::power_eigensolution(negative, Q, L);
    REQUIRE(count == 2);
    check_eigenpairs(negative, Q, L, count, 1e-5);
    const auto negative_values = stored_eigenvalues(L, count);
    CHECK(negative_values[0] == doctest::Approx(-2.0).epsilon(1e-5));
    CHECK(negative_values[1] == doctest::Approx(-1.0).epsilon(1e-5));

    Mat4x4d mixed(0);
    mixed[0][0] = 5;
    mixed[1][1] = -3;
    mixed[2][2] = 1;
    Mat4x4d Q4, L4;
    count = CGLA::power_eigensolution(mixed, Q4, L4);
    REQUIRE(count == 3);
    check_eigenpairs(mixed, Q4, L4, count, 1e-5);
    const auto mixed_values = stored_eigenvalues(L4, count);
    CHECK(mixed_values[0] == doctest::Approx(-3.0).epsilon(1e-5));
    CHECK(mixed_values[1] == doctest::Approx(1.0).epsilon(1e-5));
    CHECK(mixed_values[2] == doctest::Approx(5.0).epsilon(1e-5));

    const Mat3x3d zero(0);
    count = CGLA::power_eigensolution(zero, Q, L);
    CHECK(count == 0);
}

TEST_CASE("power_eigensolution and the 2x2 solver agree on the torus tensor")
{
    // Quad-torus inner equator. Plain power iteration stopped early and
    // reported a second eigenvalue near -2.76.
    const Mat2x2d tensor(-0.254726, -2.00256, -2.00256, 0.258012);
    double smaller = 0;
    double larger = 0;
    closed_form_2x2(tensor, smaller, larger);
    CHECK(smaller == doctest::Approx(-2.01726).epsilon(1e-4));
    CHECK(larger == doctest::Approx(2.02055).epsilon(1e-4));

    Mat2x2d Q, L;
    const int direct = CGLA::eigensolution_2x2(tensor, Q, L);
    REQUIRE(direct == 2);
    CHECK(L[0][0] == doctest::Approx(larger).epsilon(1e-12));
    CHECK(L[1][1] == doctest::Approx(smaller).epsilon(1e-12));
    check_eigenpairs(tensor, Q, L, direct, 1e-8);

    const int iterated = CGLA::power_eigensolution(tensor, Q, L);
    REQUIRE(iterated == 2);
    check_eigenpairs(tensor, Q, L, iterated, 1e-5);
    const auto values = stored_eigenvalues(L, iterated);
    CHECK(values[0] == doctest::Approx(smaller).epsilon(1e-5));
    CHECK(values[1] == doctest::Approx(larger).epsilon(1e-5));
}

TEST_CASE("eigensolution_2x2 handles a repeated eigenvalue")
{
    const Mat2x2d identity(1, 0, 0, 1);
    Mat2x2d Q, L;
    REQUIRE(CGLA::eigensolution_2x2(identity, Q, L) == 2);
    CHECK(L[0][0] == doctest::Approx(1.0));
    CHECK(L[1][1] == doctest::Approx(1.0));
    check_eigenpairs(identity, Q, L, 2, 1e-12);
}

TEST_CASE("principal curvatures of an axis-aligned paraboloid")
{
    const double radius = 0.2;
    std::vector<double> angles;
    std::vector<double> radii;
    for (int i = 0; i < 8; ++i) {
        angles.push_back(i * std::numbers::pi / 4.0);
        radii.push_back(radius);
    }

    HMesh::Manifold mesh;
    const Vec3d normal(0, 0, 1);
    const auto center = make_paraboloid(
        mesh, angles, radii, -0.5, 0.0, -2.0, Vec3d(0), Vec3d(1, 0, 0), Vec3d(0, 1, 0), normal);
    REQUIRE(center != HMesh::InvalidVertexID);

    const auto pc = HMesh::principal_curvatures(mesh, center);
    check_principal(pc, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 1e-6, 0.999);
    CHECK(std::abs(dot(pc.min_curv_direction, normal)) < 1e-6);
    CHECK(std::abs(dot(pc.max_curv_direction, normal)) < 1e-6);
}

TEST_CASE("principal curvatures follow the tangent frame")
{
    const Vec3d normal = normalize(Vec3d(0.2, -0.3, 1.0));
    const Vec3d axis_x = normalize(cross(Vec3d(0, 1, 0), normal));
    const Vec3d axis_y = cross(normal, axis_x);

    std::vector<double> angles;
    std::vector<double> radii;
    for (int i = 0; i < 8; ++i) {
        angles.push_back(0.2 + i * std::numbers::pi / 4.0);
        radii.push_back(0.2);
    }

    HMesh::Manifold mesh;
    const auto center = make_paraboloid(
        mesh, angles, radii, -0.5, 0.0, -2.0, Vec3d(0.4, -0.2, 0.1), axis_x, axis_y, normal);
    REQUIRE(center != HMesh::InvalidVertexID);

    const auto pc = HMesh::principal_curvatures(mesh, center);
    check_principal(pc, 0.5, 2.0, axis_x, axis_y, 1e-6, 0.999);
}

TEST_CASE("principal curvatures keep a cross term")
{
    std::vector<double> angles;
    std::vector<double> radii;
    for (int i = 0; i < 8; ++i) {
        angles.push_back(i * std::numbers::pi / 4.0);
        radii.push_back(0.2);
    }

    HMesh::Manifold mesh;
    const auto center = make_paraboloid(
        mesh, angles, radii, -1.0, 0.3, -1.0, Vec3d(0), Vec3d(1, 0, 0), Vec3d(0, 1, 0), Vec3d(0, 0, 1));
    REQUIRE(center != HMesh::InvalidVertexID);

    const auto pc = HMesh::principal_curvatures(mesh, center);
    check_principal(pc, 0.7, 1.3, Vec3d(1, 1, 0), Vec3d(1, -1, 0), 1e-6, 0.999);
}

TEST_CASE("principal curvatures of an irregular ring")
{
    const std::vector<double> angles = {0.1, 0.8, 1.7, 2.6, 3.5, 5.2, 6.0};
    const std::vector<double> radii = {0.08, 0.15, 0.06, 0.18, 0.11, 0.14, 0.09};

    HMesh::Manifold mesh;
    const auto center = make_paraboloid(
        mesh, angles, radii, -0.5, 0.0, -2.0, Vec3d(0), Vec3d(1, 0, 0), Vec3d(0, 1, 0), Vec3d(0, 0, 1));
    REQUIRE(center != HMesh::InvalidVertexID);

    const auto pc = HMesh::principal_curvatures(mesh, center);
    check_principal(pc, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 0.05, 0.98);

    const auto one_ring = HMesh::principal_curvatures(mesh, center, false);
    CHECK(pc.min_curvature == doctest::Approx(one_ring.min_curvature).epsilon(1e-12));
    CHECK(pc.max_curvature == doctest::Approx(one_ring.max_curvature).epsilon(1e-12));
}

TEST_CASE("opposite link vertices extend the paraboloid fit")
{
    std::vector<double> angles;
    std::vector<double> radii;
    for (int i = 0; i < 8; ++i) {
        angles.push_back(i * std::numbers::pi / 4.0);
        radii.push_back(0.2);
    }

    HMesh::Manifold mesh;
    std::vector<HMesh::VertexID> outers;
    const auto center = make_paraboloid(
        mesh, angles, radii, -0.5, 0.0, -2.0, Vec3d(0), Vec3d(1, 0, 0), Vec3d(0, 1, 0), Vec3d(0, 0, 1),
        8, &outers);
    REQUIRE(center != HMesh::InvalidVertexID);
    REQUIRE(outers.size() == 8);

    const auto widened = HMesh::principal_curvatures(mesh, center, true);
    const auto one_ring = HMesh::principal_curvatures(mesh, center, false);
    check_principal(widened, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 1e-6, 0.999);
    check_principal(one_ring, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 1e-6, 0.999);

    mesh.positions[outers.front()] += Vec3d(0, 0, 0.05);
    const auto shifted = HMesh::principal_curvatures(mesh, center, true);
    const auto ignored = HMesh::principal_curvatures(mesh, center, false);
    const double near_outer = std::abs(shifted.max_curvature - 2.0);
    CHECK(near_outer > 0.01);
    check_principal(ignored, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 1e-6, 0.999);

    // The same height error, twice as far, lies several Gaussian widths out.
    // The 1-ring sets that width, so this sample must move the fit less.
    const Vec3d placed = mesh.positions[outers.front()] - Vec3d(0, 0, 0.05);
    const double x = 2.0 * placed[0];
    const double y = 2.0 * placed[1];
    mesh.positions[outers.front()] = Vec3d(x, y, -0.25 * x * x - y * y + 0.05);
    const auto distant = HMesh::principal_curvatures(mesh, center, true);
    CHECK(std::abs(distant.max_curvature - 2.0) < near_outer);

    // The same far sample, without the Gaussian weight, pulls the fit harder.
    const auto unweighted = HMesh::principal_curvatures(mesh, center, true, false);
    const double distant_err = std::max(std::abs(distant.min_curvature - 0.5),
                                         std::abs(distant.max_curvature - 2.0));
    const double unweighted_err = std::max(std::abs(unweighted.min_curvature - 0.5),
                                            std::abs(unweighted.max_curvature - 2.0));
    CHECK(unweighted_err > distant_err);
}

TEST_CASE("a missing opposite vertex is skipped")
{
    std::vector<double> angles;
    std::vector<double> radii;
    for (int i = 0; i < 8; ++i) {
        angles.push_back(i * std::numbers::pi / 4.0);
        radii.push_back(0.2);
    }

    HMesh::Manifold mesh;
    std::vector<HMesh::VertexID> outers;
    const auto center = make_paraboloid(
        mesh, angles, radii, -0.5, 0.0, -2.0, Vec3d(0), Vec3d(1, 0, 0), Vec3d(0, 1, 0), Vec3d(0, 0, 1),
        7, &outers);
    REQUIRE(center != HMesh::InvalidVertexID);
    REQUIRE(outers.size() == 7);

    const auto pc = HMesh::principal_curvatures(mesh, center);
    check_principal(pc, 0.5, 2.0, Vec3d(1, 0, 0), Vec3d(0, 1, 0), 1e-6, 0.999);
}
