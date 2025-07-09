//
// Created by Cem Akarsubasi on 5/26/25.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <nanobench.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/cleanup.h>
#include <GEL/HMesh/face_loop.h>

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
    SUBCASE("Borderline non-manifold")
    {
        // This tests a case where a mesh is not technically manifold but this case should still be 
        // supported.
        // We first create a grid of 2x2 faces 
        auto m = create_rectangular_manifold(3,3);
        // We check that this is a single connected component
        CHECK(HMesh::connected_components(m).size()==1);
        // We now remove two diagonally opposite faces
        m.remove_face(HMesh::FaceID(0));
        m.remove_face(HMesh::FaceID(3));
        // We check that there is stil a single connected component since 
        // the two remaining faces are connected by a single vertex
        CHECK(HMesh::connected_components(m).size()==1);
        // There is also just one boundary curve
        CHECK(HMesh::count_boundary_curves(m)==1);
        // Final validation.
        CHECK(validate_manifold(m));
    }
    SUBCASE("Two tetrahedra share a vertex")
    {
        // This is a simple case where two tetra share a vertex. The important thing
        // to validate is that we get two connected components and not just one
        std::vector<CGLA::Vec3d> pts = {
            CGLA::Vec3d(0,0,0), CGLA::Vec3d(1,0,0), CGLA::Vec3d(0,1,0), 
            CGLA::Vec3d(0,0,0.5), 
            CGLA::Vec3d(0,0,1), CGLA::Vec3d(1,0,1), CGLA::Vec3d(0,1,1)
        };
        auto m = HMesh::Manifold();
        m.add_face({pts[0],pts[2], pts[1]});
        m.add_face({pts[0],pts[1], pts[3]});
        m.add_face({pts[0],pts[3], pts[2]});
        m.add_face({pts[1],pts[2], pts[3]});
        m.add_face({pts[4],pts[5], pts[6]});
        m.add_face({pts[4],pts[3], pts[5]});
        m.add_face({pts[4],pts[6], pts[3]});
        m.add_face({pts[3],pts[6], pts[5]});
        HMesh::stitch_mesh(m, 1e-6);
        // THere must be two connected components
        CHECK(HMesh::connected_components(m).size()==2);
        // no boundary curves
        CHECK(HMesh::count_boundary_curves(m)==0);
        // Final validation.
        CHECK(validate_manifold(m));
    }
    SUBCASE("Non-manifold edge")
    {
        // This is a simple case where three triangles share an edge
        // this cannot be stitched, so the best to hope for is that we
        // get two connected components, i.e. two out of three triangles
        // are stitched,
        std::vector<CGLA::Vec3d> pts = {
            CGLA::Vec3d(0,0,0), 
            CGLA::Vec3d(0,0,1), 
            CGLA::Vec3d(1,0,0), 
            CGLA::Vec3d(-1,0,0),
            CGLA::Vec3d(0,1,0)
        };
        auto m = HMesh::Manifold();
        m.add_face({pts[0],pts[1], pts[2]});
        m.add_face({pts[1],pts[0], pts[3]});
        m.add_face({pts[0],pts[1], pts[4]});

        HMesh::stitch_mesh(m, 1e-6);
        // THere must be two connected components
        CHECK(HMesh::connected_components(m).size()==2);
        // no boundary curves
        CHECK(HMesh::count_boundary_curves(m)==2);
        // Final validation.
        CHECK(validate_manifold(m));
    }
    SUBCASE("Strange fin")
    {
        // We start by forming a mesh consisting of two triangles.
        std::vector<CGLA::Vec3d> pts = {
            CGLA::Vec3d(0,0,0), 
            CGLA::Vec3d(0,1,0), 
            CGLA::Vec3d(-1,0,0), 
            CGLA::Vec3d(1,0,0)
        };
        auto m = HMesh::Manifold();
        auto f0 = m.add_face({pts[0],pts[1], pts[2]});
        m.add_face({pts[1],pts[0], pts[3]});

        HMesh::stitch_mesh(m, 1e-6);
        // THere must be a single connected components
        CHECK(HMesh::connected_components(m).size()==1);
        // one boundary curve
        CHECK(HMesh::count_boundary_curves(m)==1);
        // Final validation.
        CHECK(validate_manifold(m));

        // Now extrude the shared edge.
        HMesh::HalfEdgeSet hset = {m.walker(f0).halfedge()};
        auto extruded_faces = HMesh::extrude_halfedge_set(m, hset);
        auto f = *(extruded_faces.begin());
        auto w = m.walker(f);
        size_t k = 0;
        while (extruded_faces.count(w.opp().face()) == 0) {
            w = w.next();
            ++k;
        }
        CHECK(k<4);

        // which we cannot merge since it would result in valence 1 vertices.
        bool can_we_merge = m.merge_faces(f, w.halfedge());
        CHECK(can_we_merge == false);

        // Since we have two valence two vertices, we cannot validate that the mesh
        // is fully valid since the test is too picky for that.
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