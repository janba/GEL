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
    /// Number of nearest neighbors for the cone search, clustering, the
    /// shortest-path orientation tree, and hierarchical collapse.
    /// The non-hierarchical reconstruction graph is the Delaunay graph,
    /// not this neighborhood.
    int32_t num_neighbors = 30;
    /// Delaunay radius as a multiple of the median nearest-neighbor distance.
    /// Non-hierarchical RsR uses that graph, and MST normal alignment is its
    /// tunnel-weighted spanning forest. The same multiple also drops a
    /// component edge longer than the average edge length times this value.
    double max_neighbor_dist = 20;
    /// Cross connection angle threshold in degrees. For the k-nearest
    /// neighbor graph, if the angles between the normals of two
    /// neighbors is larger than this value, the edge is culled. This
    /// can help eliminating edges resulting from noise or edges
    /// belonging to different (potentially overlapping) components.
    double max_normal_ang = 60;
    /// Minimum step threshold for handle connections. During the handle
    /// connection phase, if the shortest path between two vertices of
    /// a handle candidate is fewer than this number, the handle candidate
    /// is rejected. This can help prevent spurious local handles from
    /// being inserted.
    int32_t max_handle_dist = 50;
    /// Distance function to use. The options are the Euclidean distance
    /// which tends to work well with smoother meshes and tangent space
    /// distance which tends to work much better with noisy meshes but
    /// with substantial quality loss in smooth meshes.
    Distance dist = Distance::Euclidean;
    /// Largest boundary loop, measured in edges, that is replaced by a
    /// face after reconstruction and then ear-clipped into triangles.
    /// A non-positive value leaves every boundary loop open.
    int32_t max_hole_edges = 100;
    /// Tree orientation only. On a tough edge, a cone from the new vertex
    /// toward its parent that contains another sheet gives the new normal
    /// the opposite sign from the parent. False keeps the agree-flip on every edge.
    /// Clustering ignores this flag. Reconstruction uses the MST.
    bool cone_check = false;
};

// Algorithm

/// Convert a point cloud into a Manifold using the rotation system
/// based reconstruction method.
/// @param vertices_in vertices of the point cloud
/// @param normals_in normals of the point cloud or empty vector
/// @param opts reconstruction options
/// @return reconstructed manifold mesh
void point_cloud_to_mesh(const std::vector<CGLA::Vec3d>& vertices_in,
                         const std::vector<CGLA::Vec3d>& normals_in,
                         const RSROpts& opts,
                         HMesh::Manifold& result);

/// How estimated normals are given a consistent sign.
/// MinimumSpanningTree is what reconstruction uses. Its graph is the
/// radius-limited Delaunay graph, and its cost is the tunnel-weighted length.
/// ShortestPathTree keeps the k-nearest neighbor graph and its edge cost, and
/// the tree is a shortest path from the same kind of root. Clustering joins
/// tangent edges before gap edges. None keeps the sign from the covariance estimate.
enum class NormalOrientation {
    Clustering = 0,
    MinimumSpanningTree = 1,
    ShortestPathTree = 2,
    None = 3,
};

/// One edge of an orientation tree. `cone_invoked` is true when this edge
/// was tough enough to run the sheet cone.
struct OrientationTreeEdge {
    int parent = 0;
    int child = 0;
    bool cone_invoked = false;
};

/// One normal per input point, computed the way reconstruction does
/// when the caller supplies no normals. `opts.dist` selects the Euclidean or
/// tangent estimate. `opts.num_neighbors` sizes the shortest-path tree,
/// clustering, and the sheet cone. The minimum spanning tree uses the
/// Delaunay graph whose radius is `opts.max_neighbor_dist` times the median
/// nearest-neighbor distance.
/// `orientation` selects the minimum spanning tree, the shortest-path tree,
/// clustering, or no sign propagation. Reconstruction calls this with the
/// minimum spanning tree.
/// `opts.cone_check` selects the sheet cone on tough tree edges. Clustering
/// ignores it.
/// When `orientation_tree` is set and `orientation` is a tree, it receives
/// one edge per tree link. Clustering leaves it empty.
void point_cloud_normals(const std::vector<CGLA::Vec3d>& vertices,
                         const RSROpts& opts,
                         std::vector<CGLA::Vec3d>& normals,
                         NormalOrientation orientation = NormalOrientation::MinimumSpanningTree,
                         std::vector<OrientationTreeEdge>* orientation_tree = nullptr);

/// Reconstruct a manifold from one undirected graph.
/// The graph is split into connected components. Faces come from the
/// rotation-system checks on the edges the graph already has. `positions` are
/// written to the mesh and are also the positions used by the rotation system.
/// `normals` are unit normals, one per node.
/// @param graph dense graph, positions, and normals with matching node ids
/// @param opts reconstruction options
/// @param result reconstructed manifold mesh
void graph_to_mesh(const ReconstructionGraph& graph,
                   const RSROpts& opts,
                   HMesh::Manifold& result);

/// Convert a point cloud into a Manifold using the hierarchical collapse
/// and reexpansion method. Rotation system reconstruction is used to perform
/// the intermediate reconstruction. The collapsed points are passed to RsR,
/// which builds a radius-limited Delaunay graph.
/// @param vertices vertices of the point cloud
/// @param normals normals of the point cloud or empty vector
/// @param collapse_options collapse options
/// @param reconstruction_options reconstruction options
/// @param reexpand_options reexpansion options
/// @return reconstructed manifold mesh
void point_cloud_collapse_reexpand(
    const std::vector<CGLA::Vec3d>& vertices,
    const std::vector<CGLA::Vec3d>& normals,
    const CollapseOpts& collapse_options,
    const RSROpts& reconstruction_options,
    const ReexpandOpts& reexpand_options,
    Manifold& manifold);

/// Hierarchical collapse and reexpansion. The collapse graph is the RsR
/// neighborhood. Only edges inside the one-ring (`CollapseOpts::initial_neighbors`)
/// are contracted. The longer edges remain and are the edges RsR triangulates.
/// Normals are the averages stored during contraction, normalized to unit length.
void point_cloud_collapse_reexpand_graph(
    const std::vector<CGLA::Vec3d>& vertices,
    const std::vector<CGLA::Vec3d>& normals,
    const CollapseOpts& collapse_options,
    const RSROpts& reconstruction_options,
    const ReexpandOpts& reexpand_options,
    Manifold& manifold);

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
