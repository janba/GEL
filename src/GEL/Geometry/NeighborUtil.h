//
// Created by Cem Akarsubasi on 11/19/25.
//

/// @file NeighborUtil.h Helper functions for dealing with point cloud neighbors

#ifndef GEL_NEIGHBORUTIL_H
#define GEL_NEIGHBORUTIL_H

#include <GEL/Util/ParallelAdapters.h>

#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>

#include <GEL/CGLA/Vec.h>

namespace Geometry
{
/// @brief Calculate tangent plane projection distance
/// @param edge: the edge to be considered
/// @param this_normal: normal of one vertex
/// @param neighbor_normal: normal of another vertex
/// @return projection distance
inline double tangent_space_distance(const CGLA::Vec3d& edge, const CGLA::Vec3d& this_normal, const CGLA::Vec3d& neighbor_normal)
{
    const double euclidean_distance = edge.length();
    if (std::abs(dot(this_normal, neighbor_normal)) < std::cos(15. / 180. * M_PI))
        return euclidean_distance;
    const double neighbor_normal_length = dot(edge, neighbor_normal);
    const double normal_length = dot(edge, this_normal);
    const double projection_dist = (std::sqrt((euclidean_distance * euclidean_distance) - (normal_length * normal_length))
    + std::sqrt((euclidean_distance * euclidean_distance) - (neighbor_normal_length * neighbor_normal_length))) * 0.5;
    return projection_dist;
}

using Tree = KDTree<CGLA::Vec3d, AMGraph::NodeID>;
using Record = KDTreeRecord<CGLA::Vec3d, AMGraph::NodeID>;

struct NeighborInfo {
    AMGraph::NodeID id;
    double distance;

    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
};

//using NeighborArray = Util::InplaceVector<NeighborInfo, 192>;
using NeighborArray = std::vector<NeighborInfo>;
using NeighborMap = std::vector<NeighborArray>;

template <typename Indices>
Tree build_kd_tree_of_indices(const std::vector<CGLA::Vec3d>& vertices, const Indices& indices)
{
    Tree kd_tree;
    for (const auto idx : indices) {
        kd_tree.insert(vertices.at(idx), idx);
    }
    kd_tree.build();
    return kd_tree;
}

/// @brief k nearest neighbor search
/// @param query: the coordinate of the point to be queried
/// @param kdTree: kd-tree for knn query
/// @param num: number of nearest neighbors to be queried
/// @param neighbors: [OUT] indices of k nearest neighbors
inline void knn_search(const CGLA::Vec3d& query, const Tree& kdTree, const int num, NeighborArray& neighbors)
{
    // It might be a better idea for the caller to handle this to reduce some clutter
    std::vector<Record> records = kdTree.m_closest(num, query, INFINITY);
    std::sort_heap(records.begin(), records.end());

    for (auto record : records) {
        neighbors.emplace_back(record);
    }
}

template <typename Indices>
auto calculate_neighbors(
    Util::detail::IExecutor& pool,
    const std::vector<CGLA::Vec3d>& vertices,
    const Indices& indices,
    const Tree& kdTree,
    const int k,
    NeighborMap&& neighbors_memoized = NeighborMap())
    -> NeighborMap
{
    if (neighbors_memoized.empty()) {
        neighbors_memoized = NeighborMap(indices.size());
        for (auto& neighbors : neighbors_memoized) {
            neighbors.reserve(k + 2);
        }
    } else if (neighbors_memoized.at(0).capacity() < k) {
        for (auto& neighbors : neighbors_memoized) {
            neighbors.reserve(k + 2);
        }
    }

    auto cache_kNN_search = [&kdTree, k, &vertices](auto index, auto& neighbor) {
        auto vertex = vertices.at(index);
        knn_search(vertex, kdTree, k, neighbor);
    };
    Util::detail::Parallel::foreach2(pool, indices, neighbors_memoized, cache_kNN_search);
    return neighbors_memoized;
}

inline auto calculate_neighbors(
    Util::detail::IExecutor& pool,
    const std::vector<CGLA::Vec3d>& vertices,
    const Tree& kdTree,
    const int k,
    NeighborMap&& neighbors_memoized = NeighborMap())
    -> NeighborMap
{
    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    return calculate_neighbors(pool, vertices, indices, kdTree, k, std::move(neighbors_memoized));
}
}

#endif //GEL_NEIGHBORUTIL_H
