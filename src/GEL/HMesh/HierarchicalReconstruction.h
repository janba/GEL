//
// Created by Cem Akarsubasi on 7/11/25.
// Some Collapse data structures and functions that can be separated from RsR

#ifndef GEL_HMESH_COLLAPSE_H
#define GEL_HMESH_COLLAPSE_H

#include <GEL/HMesh/Manifold.h>

#include <vector>
#include <numbers>

namespace HMesh::RSR
{

/// Type of distance function to use during reconstruction
/// Euclidean tends to perform better on smooth meshes while
/// Tangent distance performs better on noisy meshes.
enum class Distance {
    Euclidean,
    Tangent,
};

/// Options struct for the collapse phase of hierarchical reconstruction
struct CollapseOpts {
    /// Number of collapse iterations to run. During each iteration, a percentage
    /// of vertices determined by the reduction_per_iteration will be removed from
    /// the point cloud. Removing more than about 90% of the vertices (4 iterations
    /// at 0.5 reduction) will typically provide diminishing returns in terms of
    /// performance improvements and may affect the reconstruction results negatively
    /// as the RSR step may fail to cope with the lack of input points.
    ///
    /// The separation of this from reduction_per_iteration is to allow for future
    /// spatial splitting schemes for parallelization.
    size_t max_iterations = 4;
    /// Ratio of vertices that will be _removed_ after each iteration. See the
    /// documentation for max_iterations for details.
    double reduction_per_iteration = 0.5;
    /// Number of initial neighbors to seed for the collapse graph. All collapses
    /// should happen within the one ring of a given vertex which under most
    /// circumstances should be the nearest 4-6 vertices. Using values larger than
    /// 6 will typically not provide any quality improvement and will substantially
    /// slow down the collapse phase.
    size_t initial_neighbors = 5;
    /// Maximum number of collapses to perform. If set to a value other than 0, this
    /// will stop the collapse phase after that number of collapses are performed.
    /// For debug purposes.
    size_t max_collapses = 0;
    /// Distance function to be used for smoothing and normal estimation. Not directly
    /// used by the collapse phase.
    Distance distance = Distance::Euclidean;
};

/// Debug options for the reconstruction phase.
/// Use operator| to combine debug options.
enum ReexpandDebugOpts {
    /// No debug options
    RE_NONE            = 0,
    /// Debug prints for each iteration. Very verbose.
    RE_ITERATION       = 1 << 0,
    /// Debug prints for each vertex split. Very verbose.
    RE_SPLITS          = 1 << 1,
    /// Debug prints for the chosen vertex split. Relatively verbose.
    RE_SPLIT_RESULTS   = 1 << 2,
    /// When enabled, the reconstruction will check to see if any reexpansion
    /// creates a dihedral angle that is worse than the "angle_bad_threshold". If
    /// found, information will be printed about the expansion and a somewhat
    /// prominent floating triangle will be inserted to where the bad expansion
    /// occurred.
    RE_ERRORS          = 1 << 3,
    /// Debug prints about whether an edge flip is performed during the edge
    /// crossing checks. Very verbose.
    RE_CROSSING_FLIP   = 1 << 4,
    /// When enabled, a floating triangle will be inserted for every performed vertex
    /// split. The points of the vertex are placed to both the reexpansion points and
    /// the shared point before the reexpansion. This can be used to visualize the
    /// progression of the reexpansion phase.
    RE_MARK_SPLITS     = 1 << 6,
};

struct ReexpandDebug {
    bool early_stop_at_error = false;
    ReexpandDebugOpts debug_mask = RE_NONE;
    int stop_at_error = 0;
    size_t stop_at_iteration = 0;
    double angle_bad_threshold = 1.25;
};

/// Options struct for the reexpansion phase of hierarchical reconstruction
struct ReexpandOpts {
    /// Whether to perform the reexpansion at all. You probably want to disable
    /// this if you want to get a decimated reexpansion or want to debug the
    /// collapse phase.
    bool enabled = true;
    /// Number of iterations of the optimization algorithm to run after a vertex
    /// split. Setting this to 0 will substantial reduce the reconstruction
    /// quality. Setting this to more than 1 will typically not improve the
    /// quality much but substantially increase the runtime.
    unsigned refinement_iterations = 1;
    /// During the refinement phase, if an edge flip candidate has a corner
    /// angle in radians that is _lower_ than this value, that edge will flipped
    /// (if such a flip will not degenerate the mesh). Values around 22.5 degrees
    /// seem to provide the best results. This option will most likely be removed
    /// in the future.
    double refinement_angle_threshold = std::numbers::pi / 180.0 * 22.5;
    /// How much weight should be given to the minimum angles during the vertex
    /// splitting phase. Giving some weights to minimum angles seem to avoid a
    /// small number of degenerate cases but giving large weights to this will
    /// typically not improve quality and reducing it might improve quality instead.
    /// This option will most likely be removed in the future.
    double min_angle_weight = 1;
    /// If a triangle with a minimum angle in radians lower than min_angle_threshold
    /// is created during the vertex splitting phase, that expansion is penalized with
    /// min_angle_threshold_penalty.
    double min_angle_threshold = std::numbers::pi / 180.0 * 3.0;
    /// How much to penalize vertex splits that create triangles with minimum
    /// angles below min_angle_threshold. The value is absolute so pretty much
    /// any nonzero value will rule out those triangles.
    double min_angle_threshold_penalty = 0.1;
    /// Debug options
    ReexpandDebug debug_opts = ReexpandDebug();
};

/// Information about a single collapse.
struct SingleCollapse {
    /// Old coordinates of the active point
    CGLA::Vec3d active_point_coords;
    /// Old coordinates of the latent point
    CGLA::Vec3d latent_point_coords;
    /// Current coordinates of the active point
    CGLA::Vec3d v_bar;
};

/// A point cloud consists of pairs of coordinates and their
/// associated normal vectors.
struct PointCloud {
    std::vector<CGLA::Vec3d> points;
    std::vector<CGLA::Vec3d> normals;
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

/// Perform the collapse phase of the hierarchical reconstruction. The vertices must be provided
/// with the correct associated normals.
/// @param vertices Vertices of the point cloud to collapse
/// @param normals Normals of the point cloud to collapse, must be the same length as the vertices
/// @param opts Collapse options
/// @return A pair consisting of the collapse information and the new point cloud after the collapse
auto collapse_points(
    const std::vector<CGLA::Vec3d>& vertices,
    const std::vector<CGLA::Vec3d>& normals,
    const CollapseOpts& opts = CollapseOpts()
) -> std::pair<Collapse, PointCloud>;

/// Perform the reexpansion phase of the hierarchical reconstruction, using a manifold acquired
/// from the point set based reconstruction of the collapsed point cloud and the collapse
/// information acquired during the collapse phase.
/// @param manifold The manifold to perform the reexpansion.
/// @param collapse The collapse information acquired from the collapse phase.
/// @param opts Reexpansion options
void reexpand_points(
    Manifold& manifold,
    const Collapse& collapse,
    const ReexpandOpts& opts = ReexpandOpts());
} // namespace HMesh

#endif // GEL_HMESH_COLLAPSE_H
