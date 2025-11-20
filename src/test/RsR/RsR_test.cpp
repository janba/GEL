//
// Created by Cem Akarsubasi on 4/15/25.
//

#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/HMesh.h>
#include "../common/RawObj.h"

#include <filesystem>
#include <ranges>
#include <array>
#include <string_view>

#include <nanobench.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

using HMesh::RSR::point_cloud_to_mesh;
using namespace HMesh::RSR;
using namespace CGLA;
using NodeID = size_t;

static constexpr auto FILE_CAPITAL_A = "../../../../data/PointClouds/Capital_A.obj";
// Not included:
// static constexpr auto FILE_BUNNY_SIMPLE = "../../../../data/bunny_with_normals.obj";
static constexpr auto FILE_BUNNY_SIMPLE_NO_NORMALS = "../../../../data/bunny.obj";

// Not included:
//static constexpr auto FILE_BUNNY_COMPLEX = "../../../../data/PointClouds/bun_complete.obj";
static constexpr auto FILE_THINGY = "../../../../data/thingy.obj";
static constexpr auto FILE_AS = "../../../../data/as.obj";

template <typename... Args>
constexpr auto make_array(Args... args)
{
    constexpr auto size = sizeof...(args);
    return std::array<std::string_view, size>{args...};
}

static constexpr auto QUICK_TEST_FILES = make_array(
    FILE_CAPITAL_A,
    //FILE_THINGY,
    FILE_AS);

constexpr auto IS_EUCLIDEAN = false;
constexpr auto K_PARAM = 30;
constexpr auto GENUS = -1;
constexpr auto R_PARAM = 20;
constexpr auto THETA = 60;
constexpr auto N_PARAM = 50;

auto test_options()
{
    RSROpts opts;
    opts.dist = IS_EUCLIDEAN ? Distance::Euclidean : Distance::Tangent;
    opts.k = K_PARAM;
    opts.genus = GENUS;
    opts.r = R_PARAM;
    opts.theta = THETA;
    opts.n = N_PARAM;
    return opts;
}

template <typename Collection>
auto indices_from(const Collection& collection) -> std::vector<size_t>
{
    const auto indices = [&collection] {
        std::vector<NodeID> temp(collection.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    return indices;
}

template <typename T>
auto indexed_select(const std::vector<T>& vec, const std::vector<size_t>& indices) -> std::vector<T>
{
    std::vector<T> result;
    result.reserve(indices.size());
    for (auto idx : indices) {
        result.push_back(vec.at(idx));
    }
    return result;
}

inline
auto angle_cost(const Vec3d e1, const Vec3d e2) -> double
{
    return std::abs(2.0 * dot(e1, e2) / (e1.length() * e2.length() + 1e-18) - 1.0);
}

inline
auto edge_cost(const Vec3d e1, const Vec3d e2) -> double
{
    const Vec3d e3 = e1 - e2;
    return angle_cost(e1, e2) * e3.length();
}

inline
auto triangle_cost(const Vec3d p1, const Vec3d p2, const Vec3d p3) -> double
{
    return edge_cost((p2 - p1), (p3 - p1)) + edge_cost((p1 - p2), (p3 - p2)) + edge_cost((p1 - p3), (p2 - p3));
}

auto manifold_cost(const HMesh::Manifold& manifold) -> double
{
    double total_cost = 0.0;
    for (const auto f : manifold.faces()) {
        std::array<HMesh::VertexID, 3> arr;
        int idx = 0;
        HMesh::circulate_face_ccw(manifold, f, [&idx, &arr](const HMesh::VertexID v) {
            if (idx < 3)
                arr[idx++] = v;
        });
        if (idx == 3) {
            const auto p1 = manifold.positions[arr[0]];
            const auto p2 = manifold.positions[arr[1]];
            const auto p3 = manifold.positions[arr[2]];
            const auto cost = triangle_cost(p1, p2, p3);
            total_cost += cost;
        }
    }
    return total_cost;
}

auto manifold_average_edge_length(const HMesh::Manifold& manifold) -> double
{
    double total = 0.0;
    for (const auto halfedge : manifold.halfedges()) {
        HMesh::Walker walker = manifold.walker(halfedge);
        auto v1 = walker.vertex();
        auto v2 = walker.opp().vertex();
        GEL_ASSERT_NEQ(v1, v2);
        auto p1 = manifold.positions[v1];
        auto p2 = manifold.positions[v2];
        total += std::sqrt(dot((p2 - p1), (p2 - p1)));
    }
    return total / static_cast<double>(manifold.no_halfedges());
}

auto manifold_cost_normalized(const HMesh::Manifold& manifold) -> double
{
    const auto number_of_faces = manifold.no_faces();
    const auto total_cost = manifold_cost(manifold);

    return total_cost / static_cast<double>(number_of_faces) / manifold_average_edge_length(manifold);
}

auto manifold_is_triangular(const HMesh::Manifold& manifold) -> bool
{
    for (const auto f : manifold.faces()) {
        int idx = 0;
        HMesh::circulate_face_ccw(manifold, f, [&idx](const HMesh::VertexID v) {
            idx++;
        });
        if (idx != 3)
            return false;
    }
    return true;
}

auto manifold_face_vert_ratio(const HMesh::Manifold& manifold) -> double
{
    return static_cast<double>(manifold.no_faces()) / static_cast<double>(manifold.no_vertices());
}

/// Estimate quality of a Manifold based on dihedral angles.
/// @details The dihedral angle between each face is used to determine how expensive it would be to "cross" between
/// those two faces. The average of this for every half-edge is used to get the total cost. Boundary edges are considered
/// to be the most expensive to cross.
///
/// This gives a quality metric in terms of how "smooth" a Manifold is. In an infinite flat plane, this function is 0.
/// However, keep in mind that this is not a generalized quality measure. The Manifold might be intentionally bumpy which
/// will return a high cost.
auto manifold_dihedral_cost(const HMesh::Manifold& manifold) -> std::pair<double, double>
{
    double total_cost = 0.0;
    double edge_scaled_cost = 0.0;
    size_t bad_edges = 0;
    size_t boundary_edges = 0;
    for (const auto half_edge : manifold.halfedges()) {
        const auto f1 = manifold.walker(half_edge).face();
        const auto f2 = manifold.walker(half_edge).opp().face();
        if (f1 != HMesh::InvalidFaceID && f2 != HMesh::InvalidFaceID) {
            const auto n1 = manifold.normal(f1);
            const auto n2 = manifold.normal(f2);
            const auto cost = std::abs(CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0);
            if (!std::isnan(cost)) {
                total_cost += cost;
                edge_scaled_cost += (cost * manifold.length(half_edge));
            }
            else
                bad_edges++;
        } else {
            // boundary edges are considered to be the most expensive
            boundary_edges ++;
            total_cost += 2.0;
        }
    }
    if (bad_edges > 0) {
        std::cerr << "dihedral cost found " << bad_edges << " bad edges\n";
    }
    std::cout << "dihedral cost found " << boundary_edges << " boundary edges\n";
    return std::make_pair(total_cost / static_cast<double>(manifold.no_halfedges()), edge_scaled_cost / static_cast<double>(manifold.no_halfedges()));
}

struct ValencyHistogram {
    static constexpr size_t max_tracked_valency = 9;
    size_t total = 0;
    std::array<size_t, max_tracked_valency> common_valency{0};
    size_t high = 0;

    friend std::ostream& operator<<(std::ostream& os, const ValencyHistogram& hist);
};

std::ostream& operator<<(std::ostream& os, const ValencyHistogram& hist)
{
    os << "Total: " << hist.total << "\n";
    constexpr auto draw_factor = 40;
    for (size_t valency = 2; valency < hist.common_valency.size(); ++valency) {
        const auto ratio = static_cast<double>(hist.common_valency[valency]) / static_cast<double>(hist.total);
        const auto blocks = static_cast<int>(ratio * draw_factor);
        os << valency << ":  ";
        for (int i = 0; i < blocks; ++i) {
            os << "#";
        }
        for (int i = blocks; i < draw_factor; ++i) {
            os << " ";
        }
        os << " (" << (ratio * 100.0) << "%)\n";
    }
    const auto ratio = static_cast<double>(hist.high) / static_cast<double>(hist.total);
    const auto blocks = static_cast<int>(ratio * draw_factor);
    os << hist.common_valency.size() << "+: ";
    for (int i = 0; i < blocks; ++i) {
        os << "#";
    }
    for (int i = blocks; i < draw_factor; ++i) {
        os << " ";
    }
    os << " (" << (ratio * 100.0) << "%)\n";
    return os;
}

auto manifold_valency_histogram(const HMesh::Manifold& manifold) -> ValencyHistogram
{
    ValencyHistogram hist;
    for (const auto v : manifold.vertices()) {
        if (const auto valency = manifold.valency(v); valency > 8) {
            hist.high++;
        } else {
            hist.common_valency[valency]++;
        }
    }
    hist.total = manifold.no_vertices();
    return hist;
}

auto manifold_valency_histogram_non_boundary(const HMesh::Manifold& manifold) -> ValencyHistogram
{
    ValencyHistogram hist;
    size_t non_boundary_vertices = 0;
    for (const auto v : manifold.vertices()) {
        if (!manifold.boundary(v)) {
            non_boundary_vertices++;
            if (const auto valency = manifold.valency(v); valency > 8) {
                hist.high++;
            } else {
                hist.common_valency[valency]++;
            }
        }
    }
    hist.total = non_boundary_vertices;
    return hist;
}

// Test functions begin


auto test_reconstruct_new(const std::string_view file_name, const RSROpts& opts) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
        << "Begin new function\n";
    auto input_maybe = Util::read_raw_obj(file_name);
    if (!input_maybe) {
        WARN_FALSE("File not found.");
        return std::nullopt;
    }
    const auto input = *std::move(input_maybe);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output = point_cloud_to_mesh(input.vertices, {}, opts);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

auto point_cloud_to_mesh_legacy(std::vector<CGLA::Vec3d> const& points, const std::vector<CGLA::Vec3d>& normals, const RSROpts& opts) -> HMesh::Manifold
{
    auto points_copy = points;
    auto normals_copy = normals;
    HMesh::Manifold output;
    reconstruct_single(output, points_copy, normals_copy, opts.dist == Distance::Euclidean, opts.genus, opts.k,
                       opts.r, opts.theta, opts.n);
    return output;
}


auto test_reconstruct_legacy(const std::string_view file_name, const RSROpts& opts) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
        << "Begin original function\n";
    auto input_maybe = Util::read_raw_obj(file_name);
    if (!input_maybe) {
        WARN_FALSE("File not found.");
        return std::nullopt;
    }
    auto input = Util::to_triangle_mesh(*input_maybe);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output;
    reconstruct_single(output,
        input.vertices,
        input.normals,
        opts.dist == Distance::Euclidean, opts.genus, opts.k,
                       opts.r, opts.theta, opts.n);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

auto test_reconstruct_collapse_reexpand(const std::string_view file_name, const CollapseOpts& collapse_opts,
                                        const RSROpts& rsr_opts, const ReexpandOpts& reexpand) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
        << "Begin new function\n";
    auto input_maybe = Util::read_raw_obj(file_name);
    if (!input_maybe) {
        WARN_FALSE("File not found.");
        return std::nullopt;
    }
    const auto input = Util::to_triangle_mesh(*input_maybe);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output = point_cloud_collapse_reexpand(input.vertices, input.normals, collapse_opts, rsr_opts,
                                                           reexpand);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

auto reconstruct_assertions(const HMesh::Manifold& manifold) -> void
{
    // Manifold properties:
    // Triangular
    // Approx. 2 face per vertex (if good reconstruction)
    // Each vertex around 6 edges
    CHECK_EQ(manifold_is_triangular(manifold), true);
    auto ratio = manifold_face_vert_ratio(manifold);
    CHECK_GT(ratio, 1.5);
    CHECK_LT(ratio, 2.25);
    auto faces = manifold.no_faces();
    auto cost = manifold_cost(manifold);
    auto average_edge = manifold_average_edge_length(manifold);
    auto normed_cost = cost / static_cast<double>(faces) / average_edge;
    auto valencies = manifold_valency_histogram(manifold);
    auto valencies_nb = manifold_valency_histogram_non_boundary(manifold);
    auto [dihedral_cost1, dihedral_cost2] = manifold_dihedral_cost(manifold);
    std::cout << "ratio: " << ratio << std::endl;
    std::cout << "total cost: " << cost << std::endl;
    std::cout << "num faces : " << faces << std::endl;
    std::cout << "average edge length: " << average_edge << std::endl;
    std::cout << "normed cost: " << normed_cost << std::endl;
    std::cout << "Valency histogram (all vertices): " << std::endl;
    std::cout << valencies;
    std::cout << "Valency histogram (non-boundary vertices): " << std::endl;
    std::cout << valencies_nb;
    std::cout << "Dihedral cost raw: " << dihedral_cost1 << std::endl;
    std::cout << "Dihedral cost adj: " << dihedral_cost2 << std::endl;
}

template <typename Func>
void test_reconstruct(Func&& f, const bool save, const bool all = false)
    requires std::is_same_v<decltype(f(std::declval<std::string_view>(), std::declval<const HMesh::RSR::RSROpts&>())),
                            std::optional<HMesh::Manifold>>
{
    const auto test_files = std::ranges::subrange(QUICK_TEST_FILES);

    const auto opts = test_options();
    for (const auto file : test_files) {
        std::filesystem::path p = file;
        auto opts_neighbors = opts;
        opts_neighbors.dist = Distance::Tangent;

        auto case_name = p.stem().concat("_neighbors");
        SUBCASE(case_name.c_str()) {
            std::optional<HMesh::Manifold> manifold = f(file, opts_neighbors);
            if (manifold.has_value()) {
                auto out_path = p.stem().concat("_neighbors").concat(".obj");
                if (save)
                    HMesh::obj_save(out_path, *manifold);
                reconstruct_assertions(*manifold);
            }
        }
    }
    for (const auto file : test_files) {
        std::filesystem::path p = file;

        auto opts_euclidean = opts;
        opts_euclidean.dist = Distance::Euclidean;

        auto case_name = p.stem().concat("_euclidean");
        SUBCASE(case_name.c_str()) {
            std::optional<HMesh::Manifold> manifold = f(file, opts_euclidean);
            if (manifold.has_value()) {
                auto out_path = p.stem().concat("_euclidean").concat(".obj");
                if (save)
                    HMesh::obj_save(out_path, *manifold);
                reconstruct_assertions(*manifold);
            }
        }
    }
}

TEST_CASE("reconstruct legacy")
{
    test_reconstruct(test_reconstruct_legacy, false);
}

TEST_CASE("reconstruct")
{
    test_reconstruct(test_reconstruct_new, true);
}

TEST_CASE("reconstruct_collapse_reexpand")
{
    auto l = []<typename T0, typename T1>(T0&& PH1, T1&& PH2) {
        return test_reconstruct_collapse_reexpand(std::forward<T0>(PH1), CollapseOpts(), std::forward<T1>(PH2), ReexpandOpts());
    };
    test_reconstruct(l, true, true);
}

// TEST_CASE("reconstruct_debug")
// {
//     auto collapse_opts = CollapseOpts{};
//     auto opts_euclidean = test_options();
//     auto file = FILE_BUNNY_SIMPLE_NO_NORMALS;
//     collapse_opts.max_iterations = 3;
//     collapse_opts.initial_neighbors = 8;
//     collapse_opts.reduction_per_iteration = 0.50;
//     collapse_opts.max_collapses = 0;
//
//     auto reexpand_opts = ReexpandOptions();
//     reexpand_opts.enabled = true;
//     reexpand_opts.debug_print = false;
//     reexpand_opts.angle_threshold = std::numbers::pi / 180.0 * 9;
//     reexpand_opts.angle_threshold_penalty = 0.1;
//     reexpand_opts.angle_factor *= 1;
//     reexpand_opts.early_stop_at_error = false;
//     reexpand_opts.stop_at_error = 0;
//     reexpand_opts.debug_mask; // |= RE_ERRORS; //  | RE_FIRST_FLIP | RE_SECOND_FLIP | RE_ITERATION;
//     reexpand_opts.stop_at_iteration = 0;
//
//
//     opts_euclidean.dist = Distance::Euclidean;
//     auto manifold = test_reconstruct_collapse_reexpand(file, collapse_opts, opts_euclidean, reexpand_opts);
//     CHECK(manifold.has_value());
//     if (manifold.has_value()) {
//         HMesh::obj_save("debug_obj_reexpand_euclidean.obj", *manifold);
//         reconstruct_assertions(*manifold);
//     }
//     // auto opts_neighbors = test_options();
//     // opts_neighbors.dist = Distance::Tangent;
//     // manifold = test_reconstruct_collapse_reexpand(file, collapse_opts, opts_neighbors, reexpand_opts);
//     // CHECK(manifold.has_value());
//     // if (manifold.has_value()) {
//     //     HMesh::obj_save("debug_obj_reexpand_neighbors.obj", *manifold);
//     //     reconstruct_assertions(*manifold);
//     // }
// }