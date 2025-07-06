//
// Created by Cem Akarsubasi on 5/26/25.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <nanobench.h>

#include <GEL/HMesh/Manifold.h>

// Quick and dirty validation
bool validate_manifold(const HMesh::Manifold &m) {
    HMesh::VertexSet vs;
    HMesh::HalfEdgeSet hs;
    HMesh::FaceSet fs;
    return find_invalid_entities(m, vs, hs, fs);
}

HMesh::Manifold create_rectangular_manifold(const size_t x_size, const size_t y_size)
{
    std::vector<size_t> faces;
    std::vector<CGLA::Vec3d> vertices;
    faces.reserve((x_size - 1) * (y_size -1) * 4);
    vertices.reserve(x_size * y_size);

    for (size_t i = 0; i < y_size; ++i) {
        for (size_t j = 0; j < x_size; ++j) {
            auto x = 10.0 * static_cast<double>(j) / static_cast<double>(x_size);
            auto y = 10.0 * static_cast<double>(i) / static_cast<double>(y_size);
            vertices.emplace_back(x, y, 0.0);
        }
    }
    for (size_t i = 0; i < y_size - 1; ++i) {
        for (size_t j = 0; j < x_size - 1; ++j) {
            auto top_left = i * y_size + j;
            auto top_right = top_left + 1;
            auto bottom_left = top_left + x_size;
            auto bottom_right = bottom_left + 1;
            faces.emplace_back(top_left);
            faces.emplace_back(top_right);
            faces.emplace_back(bottom_right);
            faces.emplace_back(bottom_left);
        }
    }

    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, 4);

    return m;
}

TEST_CASE("Benchmark manifold")
{
    SUBCASE("10 x 10")
    {
        const auto x_size = 10UL;
        const auto y_size = 10UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("100 vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(validate_manifold(m));
        });
    }
    SUBCASE("100 x 100")
    {
        const auto x_size = 100UL;
        const auto y_size = 100UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("10 thousand vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(validate_manifold(m));
        });
    }
    SUBCASE("1000 x 1000")
    {
        // Approximately 60 MiB for this input
        const auto x_size = 1000UL;
        const auto y_size = 1000UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("1 million vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(validate_manifold(m));
        });
    }
    // // If Manifold::cleanup can be improved to use constant space, this should be able to run on the CI
    // // This would be the practical maximum for the largest inputs people would pass
    // SUBCASE("10000x x 10000")
    // {
    //     // Approximately 600 MiB for this input
    //     const auto x_size = 10'000UL;
    //     const auto y_size = 10'000UL;
    //     ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("100 million vertex Manifold", [&]() {
    //     auto m = create_rectangular_manifold(x_size, y_size);
    //     CHECK(validate_manifold(m));
    //     });
    // }
}