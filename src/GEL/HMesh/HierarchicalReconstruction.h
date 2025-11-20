//
// Created by Cem Akarsubasi on 7/11/25.
// Some Collapse data structures and functions that can be separated from RsR

#ifndef GEL_HMESH_COLLAPSE_H
#define GEL_HMESH_COLLAPSE_H

#include <GEL/HMesh/Manifold.h>

#include <vector>

namespace HMesh::RSR
{
enum class Distance {
    Euclidean,
    Tangent,
};

struct CollapseOpts {
    size_t max_iterations = 4;
    double reduction_per_iteration = 0.5;
    size_t initial_neighbors = 5;
    size_t max_collapses = 0;
    Distance distance = Distance::Euclidean;
};

enum DebugMask {
    RE_NONE            = 0,
    RE_ITERATION       = 1 << 0,
    RE_SPLITS          = 1 << 1,
    RE_SPLIT_RESULTS   = 1 << 2,
    RE_ERRORS          = 1 << 3,
    RE_FIRST_FLIP      = 1 << 4,
    RE_SECOND_FLIP     = 1 << 5,
    RE_MARK_SPLITS     = 1 << 6,
    RE_MARK_BAD_SPLITS = 1 << 7,
};

struct ReexpandOpts {
    bool enabled = true;
    // DEBUG OPTIONS

    bool debug_print = false;
    bool early_stop_at_error = false;
    int debug_mask = RE_NONE;
    int stop_at_error = 0;
    size_t stop_at_iteration = 0;
    double angle_stop_threshold = 1.25;
    unsigned refinement_iterations = 1;

    double angle_factor = 1;
    /// minimum angle to allow a triangle
    double angle_threshold = std::numbers::pi / 180.0 * 3.0;

    double angle_threshold_penalty = 0.1;
    double refinement_angle_threshold = 22.5;
};

// Fat 72 bytes
struct SingleCollapse {
    /// Old coordinates of the active point
    CGLA::Vec3d active_point_coords;
    /// Old coordinates of the latent point
    CGLA::Vec3d latent_point_coords;
    /// Current coordinates of the active point
    CGLA::Vec3d v_bar;
};

// FIXME: move this to CGLA
constexpr auto lerp(const CGLA::Vec3d& v1, const CGLA::Vec3d& v2, double t) -> CGLA::Vec3d
{
    return v1 * (1.0 - t) + v2 * t;
}

struct PointCloud {
    std::vector<CGLA::Vec3d> points;
    std::vector<CGLA::Vec3d> normals;
    //std::vector<size_t> indices;
};


/// Contains data needed for a reexpansion
struct Collapse {
private:
    std::vector<std::vector<SingleCollapse>> collapses;

    // Private constructors can only be called from friend functions
    Collapse() = default;
    explicit Collapse(std::vector<std::vector<SingleCollapse>>&& _collapses) : collapses(_collapses) {}

public:
    friend auto collapse_points(
        const std::vector<CGLA::Vec3d>& vertices,
        const std::vector<CGLA::Vec3d>& normals,
        const CollapseOpts& opts
    ) -> std::pair<Collapse, PointCloud>;

    friend void reexpand_points(
        Manifold& manifold,
        const Collapse& collapse,
        const ReexpandOpts& opts);
};

auto collapse_points(
    const std::vector<CGLA::Vec3d>& vertices,
    const std::vector<CGLA::Vec3d>& normals,
    const CollapseOpts& opts = CollapseOpts()
) -> std::pair<Collapse, PointCloud>;

void reexpand_points(
    Manifold& manifold,
    const Collapse& collapse,
    const ReexpandOpts& opts = ReexpandOpts());
} // namespace HMesh

#endif // GEL_HMESH_COLLAPSE_H
