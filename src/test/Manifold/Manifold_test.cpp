//
// Created by Cem Akarsubasi on 5/26/25.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>

#include <doctest.h>

#include <nanobench.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/cleanup.h>
#include <GEL/HMesh/face_loop.h>

// Validators

struct PointHash {
    size_t operator()(const CGLA::Vec3d& point) const
    {
        const auto h1 = std::hash<double>{}(point[0]);
        const auto h2 = std::hash<double>{}(point[1]);
        const auto h3 = std::hash<double>{}(point[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// Helper functions

HMesh::VertexID find_vertex_id(const HMesh::Manifold& m, const CGLA::Vec3d& p)
{
    for (size_t i = 0; i < m.no_vertices(); ++i) {
        if (m.positions[HMesh::VertexID(i)] == CGLA::Vec3d(0.0, 0.0, 0.0))
            return HMesh::VertexID(i);
    }
    return HMesh::InvalidVertexID;
}

CGLA::Vec3d half_edge_direction(const HMesh::Manifold& m, HMesh::HalfEdgeID h)
{
    const auto w = m.walker(h);
    const auto current = w.vertex();
    const auto opposing = w.opp().vertex();
    return CGLA::normalize(m.positions[opposing] - m.positions[current]);
}

std::array<HMesh::HalfEdgeID, 2> find_perpendicular_edges(const HMesh::Manifold& m, const HMesh::VertexID center_idx, const CGLA::Vec3d& to_insert_position)
{
    std::array edges = {HMesh::InvalidHalfEdgeID, HMesh::InvalidHalfEdgeID};
    const auto he_normal = CGLA::normalize(to_insert_position - m.positions[center_idx]);
    int i = 0;
    HMesh::circulate_vertex_ccw(m, center_idx, [&](HMesh::HalfEdgeID h) {
        const auto dir = half_edge_direction(m, h);
        const auto perpendicular = CGLA::dot(he_normal, dir);

        // in the perpendicular case, which is true for two edges in the test case, the above variable will be zero
        if (perpendicular == 0.0) { edges.at(i++) = h; }
    });
    return edges;
}

// HMesh construction functions

/// Makes the following shape:
///     1  -  3
///    /  \   / \
///   0  -  2 -  4
///    \  /  \ /
///     5  -  6
///
/// @return
HMesh::Manifold make_hex()
{
    std::vector<CGLA::Vec3d> vertices = {
        {-2.0, 0.0, 0.0},
        {-1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {2.0, 0.0, 0.0},
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0} };
    std::vector<size_t> faces = {
        2, 0, 1,
        2, 1, 3,
        2, 3, 4,
        2, 5, 0,
        2, 6, 5,
        2, 4, 6,
    };
    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, 3);
    CHECK_EQ(m.no_faces(), 6);
    CHECK_EQ(m.no_vertices(), 7);
    CHECK_EQ(m.no_halfedges(), 12 * 2);
    return m;
}

/// Makes the following shape:
///   0  -  2 -  4
///    \  /  \ /
///     5  -  6
///
/// @return
HMesh::Manifold make_half_hex()
{
    std::vector<CGLA::Vec3d> vertices = {
        {-2.0, 0.0, 0.0},
        {-1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {2.0, 0.0, 0.0},
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0} };
    std::vector<size_t> faces = {
        2, 5, 0,
        2, 6, 5,
        2, 4, 6,
    };
    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, 3);
    CHECK_EQ(m.no_faces(), 3);
    CHECK_EQ(m.no_vertices(), 5);
    CHECK_EQ(m.no_halfedges(), 7 * 2);
    return m;
}

/// Makes the following shape:
///     1  -  3
///       \   / \
///   0  -  2 -  4
///    \  /  \ /
///     5  -  6
///
/// @return
HMesh::Manifold make_crooked_hex_left()
{
    std::vector<CGLA::Vec3d> vertices = {
        {-2.0, 0.0, 0.0},
        {-1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {2.0, 0.0, 0.0},
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0} };
    std::vector<size_t> faces = {
        2, 1, 3,
        2, 3, 4,
        2, 5, 0,
        2, 6, 5,
        2, 4, 6,
    };
    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, 3);
    CHECK_EQ(m.no_faces(), 5);
    CHECK_EQ(m.no_vertices(), 7);
    CHECK_EQ(m.no_halfedges(), 11 * 2);
    return m;
}

/// Makes the following shape:
///     1  -  3
///    /  \   /
///   0  -  2 -  4
///    \  /  \ /
///     5  -  6
///
/// @return
HMesh::Manifold make_crooked_hex_right()
{
    std::vector<CGLA::Vec3d> vertices = {
        {-2.0, 0.0, 0.0},
        {-1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {2.0, 0.0, 0.0},
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0} };
    std::vector<size_t> faces = {
        2, 0, 1,
        2, 1, 3,
        2, 5, 0,
        2, 6, 5,
        2, 4, 6,
    };
    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, 3);
    CHECK_EQ(m.no_faces(), 5);
    CHECK_EQ(m.no_vertices(), 7);
    CHECK_EQ(m.no_halfedges(), 11 * 2);
    return m;
}

HMesh::Manifold make_final()
{
    std::vector<CGLA::Vec3d> vertices = {
        {-2.0, 0.0, 0.0},
        {-1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {2.0, 0.0, 0.0},
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0},
        {0.0, 0.5, 0.0}};
    std::vector<size_t> faces = {
        7, 0, 1,
        7, 1, 3,
        7, 3, 4,
        2, 5, 0,
        2, 6, 5,
        2, 4, 6,
        2, 0, 7, 4
    };
    std::vector<int> index_numbers = {
        3, 3, 3, 3, 3, 3, 4
    };
    HMesh::Manifold m;
    HMesh::build_manifold(m, vertices, faces, index_numbers);
    CHECK_EQ(m.no_faces(), 7);
    CHECK_EQ(m.no_vertices(), 8);
    //CHECK_EQ(m.no_halfedges(), 12 * 2);
    return m;
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

// Test cases

TEST_SUITE("slit_one_ring")
{
    TEST_CASE("center")
    {
        auto m = make_hex();
        const auto m_copy = m;
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto to_insert_position = CGLA::Vec3d(0.0, 0.5, 0.0);

        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);

        HMesh::valid(m);
        const auto vn = m.slit_one_ring(m.walker(edges[1]).opp().halfedge(), edges[0]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);

        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 6);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 14 * 2);
    }
}

TEST_SUITE("split_vertex")
{
    TEST_CASE("inner")
    {
        auto m = make_hex();
        const auto to_insert_position = CGLA::Vec3d(0.0, 0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);

        const auto vn = m.split_vertex(m.walker(edges[1]).opp().halfedge(), edges[0]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);

        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 8);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 15 * 2);

        HMesh::obj_save("hex_obj.obj", m);
    }

    TEST_CASE("boundary inner")
    {
        auto m = make_half_hex();
        const auto to_insert_position = CGLA::Vec3d(0.0, -0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(m.walker(edges[1]).face(), HMesh::InvalidFaceID);
        REQUIRE_NE(m.walker(edges[0]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[0]).opp().halfedge(), edges[1]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 5);
        CHECK_EQ(m.no_vertices(), 6);
        CHECK_EQ(m.no_halfedges(), 10 * 2);

        HMesh::obj_save("half_hex_obj_inner.obj", m);
    }

    TEST_CASE("boundary partial inner")
    {
        auto m = make_crooked_hex_left();
        const auto to_insert_position = CGLA::Vec3d(0.0, 0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(m.walker(edges[1]).face(), HMesh::InvalidFaceID);
        REQUIRE_NE(m.walker(edges[0]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[1]).opp().halfedge(), edges[0]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 7);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 14 * 2);

        HMesh::obj_save("crooked_hex_obj_inner.obj", m);
    }

    TEST_CASE("boundary partial outer left")
    {
        auto m = make_crooked_hex_left();
        const auto to_insert_position = CGLA::Vec3d(0.0, -0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(m.walker(edges[1]).face(), HMesh::InvalidFaceID);
        REQUIRE_NE(m.walker(edges[0]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[0]).opp().halfedge(), edges[1]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 7);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 14 * 2);

        HMesh::obj_save("crooked_hex_obj_outer_left.obj", m);
    }

    TEST_CASE("boundary partial outer left flipped")
    {
        auto m = make_crooked_hex_left();
        const auto to_insert_position = CGLA::Vec3d(0.0, -0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(m.walker(edges[1]).face(), HMesh::InvalidFaceID);
        REQUIRE_NE(m.walker(edges[0]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[1]).opp().halfedge(), m.walker(edges[0]).halfedge());
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 7);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 14 * 2);

        HMesh::obj_save("crooked_hex_obj_outer_left_flipped.obj", m);
    }

    TEST_CASE("boundary partial outer right")
    {
        auto m = make_crooked_hex_right();
        const auto to_insert_position = CGLA::Vec3d(0.0, -0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(m.walker(edges[1]).face(), HMesh::InvalidFaceID);
        REQUIRE_NE(m.walker(edges[0]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[0]).opp().halfedge(), edges[1]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 7);
        CHECK_EQ(m.no_vertices(), 8);
        CHECK_EQ(m.no_halfedges(), 14 * 2);

        HMesh::obj_save("crooked_hex_obj_outer_right.obj", m);
    }

    TEST_CASE("boundary outer")
    {
        auto m = make_half_hex();
        const auto to_insert_position = CGLA::Vec3d(0.0, 0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);
        REQUIRE_EQ(m.walker(edges[0]).face(), HMesh::InvalidFaceID);
        REQUIRE_EQ(m.walker(edges[1]).opp().face(), HMesh::InvalidFaceID);

        const auto vn = m.split_vertex(m.walker(edges[1]).opp().halfedge(), edges[0]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));
        CHECK_EQ(m.no_faces(), 5);
        CHECK_EQ(m.no_vertices(), 6);
        CHECK_EQ(m.no_halfedges(), 10 * 2);

        CHECK(HMesh::obj_save("half_hex_obj_outer.obj", m));
    }

    TEST_CASE("split reverse boundary")
    {
        auto m = make_half_hex();
        const auto to_insert_position = CGLA::Vec3d(0.0, -0.5, 0.0);
        const auto center_idx = find_vertex_id(m, CGLA::Vec3d(0.0, 0.0, 0.0));
        REQUIRE_NE(center_idx, HMesh::InvalidVertexID);
        const auto edges = find_perpendicular_edges(m, center_idx, to_insert_position);
        REQUIRE_NE(edges[0], HMesh::InvalidHalfEdgeID);
        REQUIRE_NE(edges[1], HMesh::InvalidHalfEdgeID);

        const auto vn = m.split_vertex(m.walker(edges[0]).opp().halfedge(), edges[1]);
        REQUIRE_NE(vn, HMesh::InvalidVertexID);
        m.positions[vn] = to_insert_position;
        CHECK(HMesh::valid(m));

        HMesh::obj_save("half_hex_obj_cursed.obj", m);
    }
}

TEST_SUITE("Misc.")
{
    TEST_CASE("Borderline non-manifold")
    {
        // This tests a case where a mesh is not technically manifold but this case should still be
        // supported.
        // We first create a grid of 2x2 faces
        auto m = create_rectangular_manifold(3,3);
        // We check that this is a single connected component
        CHECK_EQ(HMesh::connected_components(m).size(), 1);
        // We now remove two diagonally opposite faces
        m.remove_face(HMesh::FaceID(0));
        m.remove_face(HMesh::FaceID(3));
        // We check that there is stil a single connected component since
        // the two remaining faces are connected by a single vertex
        CHECK_EQ(HMesh::connected_components(m).size(), 1);
        // There is also just one boundary curve
        CHECK_EQ(HMesh::count_boundary_curves(m), 1);
        // Final validation.
        CHECK(HMesh::valid(m));
    }

    TEST_CASE("Two tetrahedra share a vertex")
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
        CHECK_EQ(HMesh::connected_components(m).size(), 2);
        // no boundary curves
        CHECK_EQ(HMesh::count_boundary_curves(m), 0);
        // Final validation.
        CHECK(HMesh::valid(m));
    }
    TEST_CASE("Non-manifold edge")
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
        CHECK_EQ(HMesh::connected_components(m).size(), 2);
        // no boundary curves
        CHECK_EQ(HMesh::count_boundary_curves(m), 2);
        // Final validation.
        CHECK(HMesh::valid(m));
    }
    // TEST_CASE("Strange fin")
    // {
    //     // We start by forming a mesh consisting of two triangles.
    //     std::vector<CGLA::Vec3d> pts = {
    //         CGLA::Vec3d(0,0,0),
    //         CGLA::Vec3d(0,1,0),
    //         CGLA::Vec3d(-1,0,0),
    //         CGLA::Vec3d(1,0,0)
    //     };
    //     auto m = HMesh::Manifold();
    //     auto f0 = m.add_face({pts[0],pts[1], pts[2]});
    //     auto f1 = m.add_face({pts[1],pts[0], pts[3]});

    //     HMesh::stitch_mesh(m, 1e-6);
    //     // THere must be a single connected components
    //     CHECK_EQ(HMesh::connected_components(m).size(), 1);
    //     // one boundary curve
    //     CHECK_EQ(HMesh::count_boundary_curves(m), 1);
    //     // Final validation.
    //     CHECK(HMesh::valid(m));

    //     // Now extrude the shared edge.
    //     HMesh::HalfEdgeID h0;
    //     for (auto h: m.incident_halfedges(f0)) {
    //         if (m.walker(h).opp().face() != HMesh::InvalidFaceID) {
    //             h0 = h;
    //             break;
    //         }
    //     }
    //     HMesh::HalfEdgeSet hset = {h0};
    //     auto extruded_faces = HMesh::extrude_halfedge_set(m, hset);
    //     auto f = *(extruded_faces.begin());
    //     // which we cannot merge since it would result in valence 1 vertices.
    //     bool can_we_merge = m.merge_faces(f, h0);
    //     CHECK_FALSE(can_we_merge);

    //     // Since we have two valence two vertices, we cannot validate that the mesh
    //     // is fully valid since the test is too picky for that.
    // }
}

TEST_CASE("Benchmark manifold")
{
    SUBCASE("10 x 10")
    {
        const auto x_size = 10UL;
        const auto y_size = 10UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("100 vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(HMesh::valid(m));
        });
    }
    SUBCASE("100 x 100")
    {
        const auto x_size = 100UL;
        const auto y_size = 100UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("10 thousand vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(HMesh::valid(m));
        });
    }
    SUBCASE("1000 x 1000")
    {
        // Approximately 60 MiB for this input
        const auto x_size = 1000UL;
        const auto y_size = 1000UL;
        ankerl::nanobench::Bench().unit("vertex").batch(x_size * y_size).run("1 million vertex Manifold", [&]() {
        auto m = create_rectangular_manifold(x_size, y_size);
        CHECK(HMesh::valid(m));
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
    //     CHECK(HMesh::valid(m));
    //     });
    // }
}