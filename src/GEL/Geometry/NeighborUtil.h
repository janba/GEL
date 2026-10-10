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

#include <algorithm>
#include <cstddef>
#include <new>
#include <ranges>

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
    AMGraph::NodeID id = 0;
    /// Euclidean distance. A kD-tree record stores the squared distance.
    double distance = 0;

    NeighborInfo() = default;
    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
    NeighborInfo(AMGraph::NodeID neighbor_id, double euclidean_distance) noexcept
        : id(neighbor_id), distance(euclidean_distance)
    {}
};

/// Neighbors from one search. The first 194 live in the object. A longer
/// search moves them to the heap. Adaptive normals request 192, and the
/// searches reserve two extra slots, so that request stays inline.
class NeighborArray {
    static constexpr size_t inline_capacity = 194;
    alignas(NeighborInfo) std::byte storage[inline_capacity * sizeof(NeighborInfo)]{};
    std::vector<NeighborInfo> heap;
    NeighborInfo* active = nullptr;
    size_t used = 0;
    size_t cap = inline_capacity;
    bool spilled = false;

    NeighborInfo* inline_data()
    {
        return reinterpret_cast<NeighborInfo*>(storage);
    }
    const NeighborInfo* inline_data() const
    {
        return reinterpret_cast<const NeighborInfo*>(storage);
    }

    void destroy_inline()
    {
        for (size_t i = 0; i < used; ++i)
            active[i].~NeighborInfo();
    }

    void spill(size_t new_cap)
    {
        std::vector<NeighborInfo> next;
        next.reserve(new_cap);
        for (size_t i = 0; i < used; ++i)
            next.push_back(active[i]);
        if (!spilled)
            destroy_inline();
        heap = std::move(next);
        spilled = true;
        active = heap.data();
        used = heap.size();
        cap = heap.capacity();
    }

public:
    using iterator = NeighborInfo*;
    using const_iterator = const NeighborInfo*;

    NeighborArray()
    {
        active = inline_data();
    }

    NeighborArray(const NeighborArray& other) : NeighborArray()
    {
        reserve(other.used);
        for (size_t i = 0; i < other.used; ++i)
            emplace_back(other.active[i]);
    }

    NeighborArray(NeighborArray&& other) noexcept : NeighborArray()
    {
        *this = std::move(other);
    }

    ~NeighborArray()
    {
        if (!spilled)
            destroy_inline();
    }

    NeighborArray& operator=(const NeighborArray& other)
    {
        if (this == &other)
            return *this;
        clear();
        reserve(other.used);
        for (size_t i = 0; i < other.used; ++i)
            emplace_back(other.active[i]);
        return *this;
    }

    NeighborArray& operator=(NeighborArray&& other) noexcept
    {
        if (this == &other)
            return *this;
        clear();
        if (other.spilled) {
            heap = std::move(other.heap);
            spilled = true;
            active = heap.data();
            used = heap.size();
            cap = heap.capacity();
            other.spilled = false;
            other.active = other.inline_data();
            other.used = 0;
            other.cap = inline_capacity;
        } else {
            active = inline_data();
            for (size_t i = 0; i < other.used; ++i)
                new (active + i) NeighborInfo(std::move(other.active[i]));
            used = other.used;
            other.destroy_inline();
            other.used = 0;
        }
        return *this;
    }

    [[nodiscard]] iterator begin() { return active; }
    [[nodiscard]] iterator end() { return active + used; }
    [[nodiscard]] const_iterator begin() const { return active; }
    [[nodiscard]] const_iterator end() const { return active + used; }

    [[nodiscard]] size_t size() const { return used; }
    [[nodiscard]] bool empty() const { return used == 0; }
    [[nodiscard]] size_t capacity() const { return cap; }

    NeighborInfo& operator[](size_t index) { return active[index]; }
    const NeighborInfo& operator[](size_t index) const { return active[index]; }

    NeighborInfo& back() { return active[used - 1]; }
    const NeighborInfo& back() const { return active[used - 1]; }

    void clear()
    {
        if (spilled)
            heap.clear();
        else if (used > 0)
            destroy_inline();
        spilled = false;
        active = inline_data();
        used = 0;
        cap = inline_capacity;
    }

    /// Drop entries past `n`. Callers shrink a filtered list.
    void resize(size_t n)
    {
        if (n >= used)
            return;
        if (spilled) {
            heap.resize(n);
            active = heap.data();
            used = heap.size();
            cap = heap.capacity();
        } else {
            for (size_t i = n; i < used; ++i)
                active[i].~NeighborInfo();
            used = n;
        }
    }

    void reserve(size_t n)
    {
        if (n <= cap)
            return;
        spill(n);
    }

    void emplace_back(const NeighborInfo& info)
    {
        if (!spilled && used == inline_capacity)
            spill(inline_capacity * 2);
        if (spilled) {
            heap.push_back(info);
            active = heap.data();
            used = heap.size();
            cap = heap.capacity();
        } else {
            new (active + used) NeighborInfo(info);
            ++used;
        }
    }
};

using NeighborMap = std::vector<NeighborArray>;

/// One k-nearest-neighbor search at one point, plus the distances read from it.
/// `nearest` is sorted nearest-first. Entry 0 is this point, at distance 0.
/// `search_neighborhoods` fills `nearest` and `farthest_distance`.
/// The other distances stay 0 until the matching `CloudNeighborhood::note_*` call.
struct PointNeighborhood {
    NeighborArray nearest;
    /// Euclidean distance of the last entry in `nearest`. 0 when `nearest` is empty.
    double farthest_distance = 0;
    /// Euclidean distance of the `one_ring_count`-th neighbor, counting this point as the first.
    double one_ring_radius = 0;
    /// Euclidean distance of the unfiltered neighbor at rank
    /// floor(`reconstruction_neighbor_count` * 2/3), clamped to the last neighbor.
    /// A recorded 0 can be this point's own distance. 0 before `note_reconstruction_cap`
    /// means the cap has not been recorded; `reconstruction_neighbor_count == 0` says which.
    double reconstruction_length_cap = 0;
};

/// The same k-nearest-neighbor search at every point. `points[i]` is vertex i.
/// `searched_count` is the k that was requested. A cloud smaller than k still
/// records that k; each list then holds every point.
/// `nearest` is shared only with a step whose requested k is `searched_count`.
/// `one_ring_radius` and `reconstruction_length_cap` are rank distances, so a
/// longer search can record a smaller rank.
struct CloudNeighborhood {
    std::vector<PointNeighborhood> points;
    int searched_count = 0;
    int one_ring_count = 0;
    int reconstruction_neighbor_count = 0;

    /// True when `nearest` is the exact `k` search of a cloud of this size.
    [[nodiscard]] bool has_exact_search(size_t point_count, int k) const
    {
        return k > 0 && searched_count == k && points.size() == point_count;
    }

    /// True when this search reaches the neighbor at `rank`, counting this point as rank 1.
    [[nodiscard]] bool has_rank(size_t point_count, int rank) const
    {
        return rank > 0 && searched_count >= rank && points.size() == point_count;
    }

    /// True when `nearest` is the exact `neighbor_count` search and
    /// `note_reconstruction_cap(neighbor_count)` has recorded the length cap.
    [[nodiscard]] bool has_reconstruction_cap(size_t point_count, int neighbor_count) const
    {
        return has_exact_search(point_count, neighbor_count)
            && reconstruction_neighbor_count == neighbor_count;
    }

    /// Record `one_ring_radius` from this search. The search must already reach `ring_count`.
    bool note_one_ring(int ring_count)
    {
        if (!has_rank(points.size(), ring_count))
            return false;
        one_ring_count = ring_count;
        for (PointNeighborhood& point : points) {
            if (point.nearest.empty()) {
                point.one_ring_radius = 0;
                continue;
            }
            const size_t rank = std::min(point.nearest.size(), static_cast<size_t>(ring_count)) - 1;
            point.one_ring_radius = point.nearest[rank].distance;
        }
        return true;
    }

    /// Record `reconstruction_length_cap` at floor(`neighbor_count` * 2/3).
    /// The index is clamped to the last returned neighbor, matching the length
    /// cap used after collapse. On a cloud larger than `neighbor_count` the
    /// index is inside the list.
    bool note_reconstruction_cap(int neighbor_count)
    {
        if (!has_rank(points.size(), neighbor_count))
            return false;
        reconstruction_neighbor_count = neighbor_count;
        const size_t cap_index = static_cast<size_t>(static_cast<double>(neighbor_count) * (2.0 / 3.0));
        for (PointNeighborhood& point : points) {
            if (point.nearest.empty()) {
                point.reconstruction_length_cap = 0;
                continue;
            }
            const size_t rank = std::min(cap_index, point.nearest.size() - 1);
            point.reconstruction_length_cap = point.nearest[rank].distance;
        }
        return true;
    }
};

template <typename Indices>
void build_kd_tree_of_indices(const std::vector<CGLA::Vec3d>& vertices, const Indices& indices, Tree& kd_tree)
{
    for (const auto idx : indices) {
        kd_tree.insert(vertices.at(idx), idx);
    }
    std::cout << "kdtree start building..." << std::endl;
    kd_tree.build();
    return;
}

/// @brief k nearest neighbor search
/// @param query: the coordinate of the point to be queried
/// @param kdTree: kd-tree for knn query
/// @param num: number of nearest neighbors to be queried
/// @param neighbors: [OUT] indices of k nearest neighbors
inline void knn_search(const CGLA::Vec3d& query, const Tree& kdTree, const int num, NeighborArray& neighbors)
{
    // Squared distances and ids. The position keys stay in the tree.
    auto records = kdTree.m_closest_ids(static_cast<unsigned>(num), query, INFINITY);
    std::sort_heap(records.begin(), records.end());

    for (const auto& record : records)
        neighbors.emplace_back(NeighborInfo(record.v, std::sqrt(record.d)));
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

/// One exact `neighbor_count` search at every vertex. `points[i].nearest` is that
/// search, and `farthest_distance` is its last Euclidean distance.
inline CloudNeighborhood search_neighborhoods(
    Util::detail::IExecutor& pool,
    const std::vector<CGLA::Vec3d>& vertices,
    const Tree& tree,
    const int neighbor_count)
{
    CloudNeighborhood cloud;
    cloud.searched_count = neighbor_count;
    cloud.points.resize(vertices.size());
    if (neighbor_count <= 0 || vertices.empty())
        return cloud;

    for (PointNeighborhood& point : cloud.points)
        point.nearest.reserve(static_cast<size_t>(neighbor_count) + 2);

    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    auto search_one = [&tree, neighbor_count, &vertices](auto index, PointNeighborhood& point) {
        knn_search(vertices[index], tree, neighbor_count, point.nearest);
        point.farthest_distance = point.nearest.empty() ? 0.0 : point.nearest.back().distance;
    };
    Util::detail::Parallel::foreach2(pool, indices, cloud.points, search_one);
    return cloud;
}
}

#endif //GEL_NEIGHBORUTIL_H
