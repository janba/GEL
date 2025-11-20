#ifndef GEL_RSR_EXPERIMENTAL_H
#define GEL_RSR_EXPERIMENTAL_H
#pragma once

#include <vector>

#include <GEL/HMesh/HierarchicalReconstruction.h>

/// @brief Rotation System Reconstruction
namespace HMesh::RSR
{
/// Options struct for point cloud reconstruction
struct RSROpts {
    /// Expected genus of the manifold, -1 to auto-detect. This value
    /// is used in the handle connection phase to determine how many
    /// handles should be inserted. Providing the genus might improve
    /// reconstruction quality but only if the reconstruction is of
    /// sufficiently good quality.
    int32_t genus = -1;
    /// Number of nearest neighbors to initialize the graphs with.
    /// Higher values can improve reconstruction quality, but they
    /// will also substantially increase reconstruction time. Values
    /// higher than 30 are likely to result in diminishing returns.
    int32_t k = 30;
    /// Threshold multiplier for removing overly long edges. When
    /// initializing k-nearest neighbor graphs, edges whose length are
    /// longer than the average length times r will be removed from
    /// the graph. This can prevent outliers and unwanted connections
    /// between different components of a manifold.
    double r = 20;
    /// Cross connection angle threshold in degrees. For the k-nearest
    /// neighbor graph, if the angles between the normals of two
    /// neighbors is larger than this value, the edge is culled. This
    /// can help eliminating edges resulting from noise or edges
    /// belonging to different (potentially overlapping) components.
    double theta = 60;
    /// Minimum step threshold for handle connections. During the handle
    /// connection phase, if the shortest path between two vertices of
    /// a handle candidate is fewer than this number, the handle candidate
    /// is rejected. This can help prevent spurious local handles from
    /// being inserted.
    int32_t n = 50;
    /// Distance function to use. The options are the Euclidean distance
    /// which tends to work well with smoother meshes and tangent space
    /// distance which tends to work much better with noisy meshes but
    /// with substantial quality loss in smooth meshes.
    Distance dist = Distance::Euclidean;
};

// Algorithm

/// Convert a point cloud into a Manifold using the rotation system
/// based reconstruction method.
/// @param vertices_in vertices of the point cloud
/// @param normals_in normals of the point cloud or empty vector
/// @param opts reconstruction options
/// @return reconstructed manifold mesh
auto point_cloud_to_mesh(const std::vector<CGLA::Vec3d>& vertices_in,
                         const std::vector<CGLA::Vec3d>& normals_in,
                         const RSROpts& opts) -> HMesh::Manifold;

/// Convert a point cloud into a Manifold using the hierarchical collapse
/// and reexpansion method. Rotation system reconstruction is used to perform
/// the intermediate reconstruction.
/// @param vertices vertices of the point cloud
/// @param normals normals of the point cloud or empty vector
/// @param collapse_options collapse options
/// @param reconstruction_options reconstruction options
/// @param reexpand_options reexpansion options
/// @return reconstructed manifold mesh
auto point_cloud_collapse_reexpand(
    const std::vector<CGLA::Vec3d>& vertices,
    const std::vector<CGLA::Vec3d>& normals,
    const CollapseOpts& collapse_options,
    const RSROpts& reconstruction_options,
    const ReexpandOpts& reexpand_options) -> HMesh::Manifold;

// TODO: probably remove these
namespace detail
{
    struct NormalEstimationResult {
        std::vector<CGLA::Vec3d> vertices;
        std::vector<CGLA::Vec3d> normals;
        std::vector<CGLA::Vec3d> smoothed_v;
    };

    auto point_cloud_normal_estimate(const std::vector<CGLA::Vec3d>& vertices,
                                     const std::vector<CGLA::Vec3d>& normals,
                                     bool is_euclidean) -> NormalEstimationResult;
}
} // namespace HMesh::RSR

#endif // GEL_HMesh_RsR2_hpp
