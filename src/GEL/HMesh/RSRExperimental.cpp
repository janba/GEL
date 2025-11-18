#include <GEL/Util/RSRTimer.h>
#include <GEL/Geometry/normal.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/Util/ParallelAdapters.h>
#include <ranges> // std::views
#include <algorithm>
//#include <GEL/Util/RangeTools.h>

#include "GEL/Util/InplaceVector.h"

namespace HMesh::RSR
{
using namespace detail;
using namespace Util::detail;

using Geometry::estimateNormal;
using namespace ::HMesh;
//using namespace Util::Ranges;
using Util::AttribVec;
using HMesh::Manifold;

using FaceType = std::array<NodeID, 3>;

struct FaceComparator {
    constexpr bool operator()(const std::pair<FaceType, float>& left,
                              const std::pair<FaceType, float>& right) const
    {
        return (left.second) > (right.second);
    }
};

using FacePriorityQueue = std::priority_queue<
    std::pair<FaceType, float>,
    std::vector<std::pair<FaceType, float>>,
    FaceComparator>;

struct TEdge {
    NodeID from;
    NodeID to;
};
struct PEdgeLength {
    TEdge edge;
    float length;
};

struct FaceConnectionKey {
    uint tree_id;
    uint to_tree_id;

    constexpr bool operator==(const FaceConnectionKey& other) const noexcept
    {
        return tree_id == other.tree_id && to_tree_id == other.to_tree_id;
    }
};

struct FaceConnectionKeyHasher {
    std::size_t operator()(const FaceConnectionKey& self) const noexcept
    {
        constexpr std::hash<size_t> hasher;
        return (hasher(self.tree_id) << 1) ^ hasher(self.to_tree_id);
    }
};

using FacePair = std::pair<int, FaceConnectionKey>;
using NeighborPair = std::pair<double, NodeID>;

constexpr bool edge_comparator(const PEdgeLength& l, const PEdgeLength& r)
{
    return l.length < r.length;
}

constexpr bool face_comparator(const FacePair& l, const FacePair& r)
{
    return l.first > r.first;
}

constexpr bool neighbor_comparator(const NeighborPair& l, const NeighborPair& r)
{
    return l.first > r.first;
}

/// @brief Calculate the reference vector for the rotation system
/// @param normal: normal direction for the target vertex
/// @return the reference vector
constexpr Vec3 calculate_ref_vec(const Vec3& normal)
{
    constexpr double eps = 1e-6;
    if (normal[2] == 1.0)
        return Vec3(0.0, 1.0, 0.0);;
    const double second = (normal[1] == 0.0) ? eps : normal[1];
    const auto ref_vec = Vec3(0.0, -normal[2] / second, 1.0);
    return CGLA::normalize(ref_vec);
}

/// @brief Calculate the radian given the reference vector
/// @param branch: vector of the outgoing edge
/// @param normal: normal of the root vertex
/// @param ref_vec: the reference vector
/// @return radian
double cal_radians_3d(const Vec3& branch, const Vec3& normal, const Vec3& ref_vec)
{
    const Vec3 proj_vec = branch - dot(normal, branch) * normal;
    const double proj_vec_len = proj_vec.length();
    if (proj_vec_len < 1e-8)
        return 0.;

    const Vec3 proj_ref = ref_vec - dot(normal, ref_vec) * normal;
    const auto value = std::clamp<double>(
        dot(proj_vec, proj_ref) / proj_vec_len /
        proj_ref.length(), -1.0, 1.0);
    // We clamped value to [-1,1], so radian cannot be NaN unless normal is 0
    double radian = std::acos(value);
    if (dot(cross(proj_vec, proj_ref), normal) > 0.0)
        radian = 2.0 * M_PI - radian;
    return radian;
}

/// @brief Calculate the radian in the rotation system
/// @param branch: vector of the outgoing edge
/// @param normal: normal of the root vertex
/// @return radian
double cal_radians_3d(const Vec3& branch, const Vec3& normal)
{
    const auto ref_vec = calculate_ref_vec(normal);
    return cal_radians_3d(branch, normal, ref_vec);
}

/// @brief Calculate projection distance
/// @param edge: the edge to be considered
/// @param this_normal: normal of one vertex
/// @param neighbor_normal: normal of another vertex
/// @return projection distance
double cal_proj_dist(const Vec3& edge, const Vec3& this_normal, const Vec3& neighbor_normal)
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

void filter_out_cross_connection(
    NeighborArray& neighbors,
    const std::vector<Vec3>& normals,
    const NodeID this_idx,
    const double cross_conn_thresh,
    const bool is_euclidean)
{
    const auto& this_normal = normals[this_idx];
    std::erase_if(neighbors, [&](const auto& neighbor_info) {
        const auto& neighbor_normal = normals[neighbor_info.id];
        const double cos_theta = dot(this_normal, neighbor_normal);
        const double cos_thresh =
            (is_euclidean) ? 0.0 : std::cos(cross_conn_thresh / 180. * M_PI);
        return cos_theta < cos_thresh;
    });
}

struct ConnectionLength {
    /// the distance of the longest connection each vertex involved
    double max_length = 0.0;
    /// the maximum length of connection before connecting handles (conservative connection)
    double pre_max_length = 0.0;
};

/// @brief initialize the graph and related information
/// @param vertices: vertices of the component
/// @param smoothed_v: smoothed vertices of the component
/// @param normals: normals of this component
/// @param neighbor_map
/// @param connection_lengths [OUT] Distance of each vertex's longest connection and the maximum length before connecting handles
/// @param k
/// @param is_euclidean
/// @return a light-weight graph with the essential connections for building MST
SimpGraph init_graph(
    const std::vector<Point>& smoothed_v,
    const std::vector<Vec3>& normals,
    const NeighborMap& neighbor_map,
    std::vector<ConnectionLength>& connection_lengths,
    const int k,
    const double cross_connection_threshold,
    const bool is_euclidean)
{
    SimpGraph dist_graph;
    AMGraph::NodeSet sets;
    dist_graph.reserve(smoothed_v.size(), k);
    //GEL_ASSERT_EQ(vertices.size(), smoothed_v.size());
    GEL_ASSERT_EQ(smoothed_v.size(), normals.size());

    for (int i = 0; i < smoothed_v.size(); i++) {
        auto node = dist_graph.add_node();
        sets.insert(node);
    }

    for (NodeID id_this = 0UL; id_this < smoothed_v.size(); ++id_this) {
        const Point& vertex = smoothed_v[id_this];
        const Vec3& this_normal = normals[id_this];
        auto neighbors = neighbor_map[id_this];

        connection_lengths[id_this].pre_max_length = neighbors[static_cast<size_t>(double(k) * (2.0 / 3.0))].distance;

        filter_out_cross_connection(neighbors, normals, id_this, cross_connection_threshold, is_euclidean);
        //GEL_ASSERT(neighbors.size() > 1);
        for (const auto& neighbor : neighbors | std::views::drop(1)) {
            const auto id_other = neighbor.id;
            //if (id_other <= id_this) continue;
            //GEL_ASSERT_NEQ(id_other, id_this, "Vertex connects back to its own");

            const Vec3&  neighbor_normal = normals[id_other];
            const Point& neighbor_pos = smoothed_v[id_other];
            const Vec3  edge = neighbor_pos - vertex;
            const auto weight =
                (is_euclidean)
                    ? (edge.length())
                    : (cal_proj_dist(edge, this_normal, neighbor_normal));
            if (weight > connection_lengths[id_this].max_length)
                connection_lengths[id_this].max_length = weight;

            if (weight > connection_lengths[id_other].max_length)
                connection_lengths[id_other].max_length = weight;

            if (weight < 1e-8) [[unlikely]]
                std::cerr << __func__ << " [ERROR] weight: " << weight << std::endl;

            // if (dist_graph.find_edge(id_this, id_other) != SimpGraph::InvalidEdgeID)
            //     continue;

            dist_graph.connect_nodes(id_this, id_other, weight);
        }
    }

    const auto components_vec = connected_components(dist_graph.inner(), sets);
    GEL_ASSERT(components_vec.size() == 1, "%d", components_vec.size());
    return dist_graph;
}

/// @brief Find the shortest path from one vertex to another in the graph
/// @param mst: the graph
/// @param start: the source vertex
/// @param target: the target vertex
/// @return number of steps in the shortest path
int find_shortest_path(const RSGraph& mst, const NodeID start, const NodeID target)
{
    std::queue<NodeID> q;
    std::vector<int> dist(mst.no_nodes(), -1); // Distance from the start to each node
    UnorderedSet<NodeID> visited;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        const auto u = q.front();
        q.pop();

        // If the target node is reached, stop early
        if (u == target)
            break;
        // If the node has already been visited, skip it
        if (visited.contains(u))
            continue;

        visited.insert(u);

        // Explore neighbors
        for (const auto& v : mst.neighbors_lazy(u)) {
            if (dist[v] == -1) {
                // If the node hasn't been visited
                dist[v] = dist[u] + 1; // Increment distance
                q.push(v);
            }
        }
    }

    return dist[target];
}

/// @brief weighted smoothing method using defined neighborhood with tangential distance weighted
/// @param pool thread pool
/// @param vertices: vertices of the point cloud
/// @param normals: normal of the point cloud
/// @param neighbors_map
/// @param smoothed_v: [OUT] vertices after smoothing
void weighted_smooth(
    IExecutor& pool,
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const NeighborMap& neighbors_map,
    std::vector<Point>& smoothed_v)
{
    auto _debug = __func__;
    auto lambda = [_debug, &normals, &vertices](const size_t idx, const Point& vertex, const NeighborArray& neighbors) {
        const Vec3 normal = normals[idx];

        double weight_sum = 0.;
        double amp_sum = 0.;
        double max_dist = 0.;

        struct LengthWeight {
            double vert_length = 0.0;
            double weight = 0.0;
        };
        InplaceVector<LengthWeight, 192> length_weights;
        const std::intptr_t limit = (neighbors.size() < 192) ? static_cast<intptr_t>(neighbors.size()) : 192;
        length_weights.reserve(limit);
        for (const auto& neighbor: neighbors | std::views::take(192)) {
            const Point neighbor_pos = vertices[neighbor.id];
            const Vec3 n2this = neighbor_pos - vertex;
            if (dot(normals[neighbor.id], normal) < std::cos(30. / 180. * M_PI)) {
                continue;
            }
            const double vertical = CGLA::dot(n2this, normal);
            const double n_dist = (neighbor_pos - vertex).length();

            const double tangential_square = n_dist * n_dist -
                vertical * vertical;
            const double tangential_dist = (tangential_square > 0.) ? std::sqrt(tangential_square) : 0.0;

            if (!std::isfinite(tangential_dist)) [[unlikely]] {
                std::cerr << n_dist << " " << vertical << std::endl;
                std::cerr << _debug << ": error" << std::endl;
            }

            const double weight = -tangential_dist;
            if (tangential_dist > max_dist)
                max_dist = tangential_dist;

            length_weights.emplace_back(vertical, weight);
            weight_sum += weight;
        }
        for (const auto& length_weight : length_weights) {
            amp_sum += length_weight.vert_length * (length_weight.weight + max_dist);
        }
        weight_sum += (max_dist * static_cast<double>(length_weights.size()));

        weight_sum = (weight_sum == 0.) ? 1. : weight_sum;

        amp_sum /= weight_sum;
        if (!std::isfinite(amp_sum)) [[unlikely]]
            std::cout << _debug << ": error" << std::endl;
        const Vec3 move = amp_sum * normal;
        return vertex + move;
    };
    Parallel::enumerate_map2(pool, vertices, neighbors_map, smoothed_v, lambda);
}

auto normalize_normals(std::vector<Vec3>& normals) -> void
{
    for (auto& normal : normals) {
        normal.normalize();
    }
}

void estimate_normal_no_normals_memoized(
    IExecutor& pool,
    const std::vector<Point>& vertices,
    const NeighborMap& neighbors,
    std::vector<Vec3>& normals)
{
    normals.clear();
    // Data type transfer & Cal diagonal size
    auto _debug = __func__;
    auto lambda = [&](const NeighborArray& neighbors_of_this) {
        // need id, distance and coords anyway

        auto neighbor_coords = neighbors_of_this | std::views::transform([&](const auto& neighbor) {
            return vertices[neighbor.id];
        });
        const Vec3 normal = estimateNormal(neighbor_coords);

        if (std::isnan(normal.length())) [[unlikely]] {
            std::cerr << _debug << ": error" << std::endl;
            std::cerr << neighbors_of_this.size() << std::endl;
        }
        return normal;
    };
    Parallel::map(pool, neighbors, normals, lambda);
}

/// @brief Calculate cos angle weight for correcting normal orientation
/// @param this_normal: normal of current vertex
/// @param neighbor_normal: normal of its neighbor vertex
/// @return angle weight calculated
double cal_angle_based_weight(const Vec3& this_normal, const Vec3& neighbor_normal)
{
    const double dot_pdt = std::abs(
        dot(this_normal, neighbor_normal) / (this_normal.length() * neighbor_normal.length()));
    const double dot_pdt_clamped = std::clamp<double>(dot_pdt, 0., 1.0);
    return 1.0 - dot_pdt_clamped;
}

/// Generic template for creating an MST
template </*std::derived_from<AMGraph>*/ typename TargetGraph,
    /*std::derived_from<AMGraph>*/ typename SourceGraph,
    std::invocable<TargetGraph&, NodeID> NodeInserter, // TargetGraph -> NodeID -> ()
    std::invocable<const SourceGraph&, NodeID, NodeID> DistanceFunc, // SourceGraph -> NodeID -> NodeID -> double
    std::invocable<TargetGraph&, NodeID, NodeID> EdgeInserterFunc> // TargetGraph -> NodeID -> NodeID -> ()
TargetGraph make_mst(
    const SourceGraph& g,
    NodeID root,
    NodeInserter&& node_inserter,
    DistanceFunc&& distance_function,
    EdgeInserterFunc&& edge_insertion)
{
    GEL_ASSERT_NEQ(root, AMGraph::InvalidNodeID);
    TargetGraph gn;
    struct QElem {
        double distance;
        NodeID from;
        NodeID to;

        bool operator>(const QElem& rhs)
        {
            return distance > rhs.distance;
        }
    };
    struct Bool {
        bool inner;
    };
    for (auto n : g.node_ids())
        node_inserter(gn, n);
    std::vector<Bool> in_tree(gn.no_nodes(), Bool(false));
    in_tree[root].inner = true;

    std::priority_queue<QElem, std::vector<QElem>, std::greater<>> queue;
    for (auto neighbor : g.neighbors_lazy(root)) {
        const auto distance = distance_function(g, neighbor, root);
        queue.emplace(distance, root, neighbor);
    }

    while (!queue.empty()) {
        auto [distance, id1, id2] = queue.top();
        queue.pop();

        if (!in_tree[id2].inner) {
            in_tree[id2].inner = true;

            edge_insertion(gn, id1, id2);

            GEL_ASSERT(std::ranges::distance(g.neighbors_lazy(id2)) > 0);
            for (auto neighbor : g.neighbors_lazy(id2)
                | std::views::filter([&](auto&& nb){return !in_tree[nb].inner;})) {
                const auto distance2 = distance_function(g, neighbor, id2);
                queue.emplace(distance2, id2, neighbor);
            }
        }
    }
    return gn;
}

RSGraph minimum_spanning_tree(
    const SimpGraph& g,
    const NodeID root,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    const bool is_euclidean)
{
    GEL_ASSERT_NEQ(root, AMGraph::InvalidNodeID);
    auto gn = make_mst<RSGraph>(
        g, root,
        [&](RSGraph& graph, NodeID n) { graph.add_node(vertices[n], normals[n]); },
        [&](const SimpGraph& graph, NodeID n, NodeID m) {
            return CGLA::sqr_length(vertices[m] - vertices[n]);
        },
        [&](RSGraph& graph, NodeID id1, NodeID id2) {
            graph.add_edge(id1, id2);
        }
    );
    // by the definition of an mst
    GEL_ASSERT_EQ(gn.m_vertices.size(), gn.m_edges.size() + 1);
    return gn;
}

SimpGraph minimum_spanning_tree(const SimpGraph& g, NodeID root)
{
    auto gn = make_mst<SimpGraph>(g, root,
                                  [](SimpGraph& gn, NodeID n) { gn.add_node(); },
                                  [](const SimpGraph& g, NodeID n, NodeID m) { return g.get_weight(n, m); },
                                  [](SimpGraph& gn, NodeID n, NodeID m) { gn.connect_nodes(n, m); });
    return gn;
}

/// @brief Determine the normal orientation
/// @param pool Thread pool
/// @param kdTree
/// @param in_smoothed_v
/// @param normals: [OUT] normal of the point cloud with orientation corrected
/// @param k
void correct_normal_orientation(
    IExecutor& pool,
    const Tree& kdTree,
    const std::vector<Point>& in_smoothed_v,
    std::vector<Vec3>& normals,
    const int k)
{
    /// The graph has the angles as weights
    const auto [g_angle, sets] = [&] {
        SimpGraph g_angle_temp;
        AMGraph::NodeSet sets_temp;
        // sets_temp.reserve(in_smoothed_v.size()); // TODO: flat_hash_set can be reserved

        for (int i = 0; i < in_smoothed_v.size(); i++) {
            sets_temp.insert(g_angle_temp.add_node());
        }

        const auto all_neighbors = calculate_neighbors(pool, in_smoothed_v, kdTree, k);

        // Init angle based graph
        for (int i = 0; i < in_smoothed_v.size(); i++) {
            const auto& neighbors = all_neighbors[i];
            const auto& this_normal = normals[i];

            for (const auto neighbor : neighbors
                    | std::views::drop(1)
                    | std::views::filter([&i, &g_angle_temp](auto&& nb)
                        {return g_angle_temp.find_edge(i, nb.id) == AMGraph::InvalidEdgeID;})
                        ) {

                //if (g_angle_temp.find_edge(i, neighbor.id) != AMGraph::InvalidEdgeID)
                //    continue;
                const auto& neighbor_normal = normals[neighbor.id];
                const double angle_weight = cal_angle_based_weight(this_normal, neighbor_normal);

                g_angle_temp.connect_nodes(i, neighbor.id, angle_weight);
            }
        }

        return std::make_tuple(g_angle_temp, sets_temp);
    }();

    const auto components_vec = connected_components(g_angle.inner(), sets);

    // The number of components and their relative sizes is inconsistent, don't parallelize this
    for (const auto& i : components_vec) {
        NodeID root = *i.begin();
        SimpGraph mst_angle = minimum_spanning_tree(g_angle, root);

        /// A trivial wrapper over bool to avoid the std::vector<bool> specialization
        struct Boolean {
            bool inner;
        };
        auto visited_vertex = std::vector(g_angle.no_nodes(), Boolean{false});

        // This uses the MST to visit every node
        // Start from the root
        std::queue<NodeID> to_visit;
        to_visit.push(root);
        while (!to_visit.empty()) {
            const NodeID node_id = to_visit.front();
            to_visit.pop();

            visited_vertex[node_id].inner = true;
            const Vec3 this_normal = normals[node_id];
            for (auto vd : mst_angle.neighbors_lazy(node_id)) {
                if (!visited_vertex[vd].inner) {
                    to_visit.push(vd);
                    const Vec3 neighbor_normal = normals[vd];
                    if (dot(this_normal, neighbor_normal) < 0) {
                        normals[vd] = -normals[vd];
                    }
                }
            }
        }
    }
}

void init_face_loop_label(RSGraph& g)
{
    constexpr NodeID start_v = 0;
    NodeID last_vertex = start_v;
    auto loop_step = 0UL;
    // TODO: unsafe access
    NodeID current_vertex = g.m_neighbors[start_v].begin()->v;
    do {
        auto& next_neighbor = g.predecessor(current_vertex, last_vertex);

        next_neighbor.tree_id = g.etf.accumulate();

        last_vertex = current_vertex;
        current_vertex = next_neighbor.v;

        loop_step++;
    } while (current_vertex != g.m_neighbors[start_v].begin()->v || last_vertex != start_v);

    std::cout << "Loop step initialization finished after " + std::to_string(loop_step) + " steps." << std::endl;
}

/// @brief Project a vector to a plane
/// @param input: Vector to be projected
/// @param normal: normal to the plane
/// @return projected Vector
Vec3 projected_vector(const Vec3& input, const Vec3& normal)
{
    const auto normal_normed = CGLA::normalize(normal);
    return input - CGLA::dot(input, normal_normed) * normal_normed;
    //return input - dot(input, normal) * normal;
}

/// @brief Check if two segments are intersecting on both planes (defined by their normal) they belong
/// @param mst: graph and vertex information
/// @param v1: 1st vertex of segment 1
/// @param v2: 2nd vertex of segment 1
/// @param v3: 1st vertex of segment 2
/// @param v4: 2nd vertex of segment 2
/// @return if they are intersecting with each other
bool is_intersecting(const RSGraph& mst, const NodeID v1, const NodeID v2, const NodeID v3, const NodeID v4)
{
    const auto& p1 = mst.m_vertices[v1].coords;
    const auto& p2 = mst.m_vertices[v2].coords;
    const auto& n1  = mst.m_vertices[v1].normal;
    const auto& n2  = mst.m_vertices[v2].normal;
    const Vec3 normal_12 = (n1 + n2) / 2.;

    const auto& p3 = mst.m_vertices[v3].coords;
    const auto& p4 = mst.m_vertices[v4].coords;
    const auto& n3  = mst.m_vertices[v3].normal;
    const auto& n4  = mst.m_vertices[v4].normal;
    const Vec3 normal_34 = (n3 + n4) / 2.;

    auto check_intersection = [](auto&& p1, auto&& p2, auto&& p3, auto&& p4, auto&& normal) -> bool {
        const Vec3 midpoint = p1 + (p2 - p1) / 2.0;
        const Vec3 edge1 = p1 - midpoint;
        const Vec3 edge2 = p3 - midpoint;
        const Vec3 edge3 = p4 - midpoint;
        const Vec3 proj_edge1 = projected_vector(edge1, normal);
        const Vec3 proj_edge2 = projected_vector(edge2, normal);
        const Vec3 proj_edge3 = projected_vector(edge3, normal);
        const Vec3 pro1 = CGLA::cross(proj_edge2, proj_edge1);
        const Vec3 pro2 = CGLA::cross(proj_edge3, proj_edge1);
        return CGLA::dot(pro1, pro2) <= 0;
    };

    // On the plane of edge 12
    if (check_intersection(p1, p2, p3, p4, normal_12) &&
        check_intersection(p3, p4, p1, p2, normal_12))
        return true;

    // On the plane of edge 34
    if (check_intersection(p1, p2, p3, p4, normal_34) &&
        check_intersection(p3, p4, p1, p2, normal_34))
        return true;

    return false;
}

/// @brief Geometry check for connection
/// @param mst: graph and vertex information
/// @param candidate: the edge to be checked
/// @param kd_tree: kd-tree for knn query
/// @return if the candidate pass the check
bool geometry_check(
    const RSGraph& mst,
    const TEdge& candidate,
    const Tree& kd_tree,
    std::vector<std::pair<Point, NodeID>>& neighbors)
{
    const NodeID& v1 = candidate.from;
    const NodeID& v2 = candidate.to;
    const Point& p1 = mst.m_vertices[v1].coords;
    const Point& p2 = mst.m_vertices[v2].coords;
    const Vec3& n1 = mst.m_vertices[v1].normal;
    const Vec3& n2 = mst.m_vertices[v2].normal;

    const Vec3 mean_normal = (n1 + n2) * 0.5;

    const Point search_center = (p1 + p2) * 0.5;
    const double radius = (p2 - p1).length() * 0.5;

    neighbors.clear();
    kd_tree.in_sphere(search_center, radius * 3.0, neighbors);

    auto in_rejection_set = [&](NodeID neighbor) -> bool {
        return
            (neighbor != v1 &&
            neighbor != v2 &&
            CGLA::dot(mst.m_vertices[neighbor].normal, mean_normal) > 0.5); // std::cos(60. / 180. * M_PI));
    };
    for (const auto neighbor : neighbors
         | std::views::values
         | std::views::filter(in_rejection_set)) {

        if (std::ranges::any_of(
            mst.m_neighbors[neighbor] |
            std::views::transform([](auto&& x) { return x.v; }) |
            std::views::filter(in_rejection_set), [&](const auto& rej_neighbor_neighbor) {
                return is_intersecting(mst, v1, v2, neighbor, rej_neighbor_neighbor);
            }
        )) return false;
    }
    return true;
}

bool geometry_check(
    const RSGraph& mst,
    const TEdge& candidate,
    const Tree& kd_tree)
{
    std::vector<std::pair<Point, NodeID>> neighbors;
    return geometry_check(mst, candidate, kd_tree, neighbors);
}

bool topology_check(const RSGraph& mst, const TEdge& candidate)
{
    const auto [this_v, neighbor] = candidate;

    // Topology check
    const auto this_v_tree = mst.predecessor(this_v, neighbor).tree_id;
    const auto neighbor_tree = mst.predecessor(neighbor, this_v).tree_id;

    if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
        return false;
    } else {
        return true;
    }
}

bool vanilla_check(
    const RSGraph& mst,
    const TEdge& candidate,
    const Tree& kd_tree,
    std::vector<std::pair<Point, NodeID>>& neighbors)
{
    const auto [this_v, neighbor] = candidate;

    // Topology check
    const auto this_v_tree = mst.predecessor(this_v, neighbor).tree_id;
    const auto neighbor_tree = mst.predecessor(neighbor, this_v).tree_id;

    if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
        return false;
    }

    return geometry_check(mst, candidate, kd_tree, neighbors);
}

bool routine_check(const RSGraph& mst, const FaceType& triangle)
{
    const Point& p1 = mst.m_vertices[triangle[0]].coords;
    const Point& p2 = mst.m_vertices[triangle[1]].coords;
    const Point& p3 = mst.m_vertices[triangle[2]].coords;
    const auto len_ui = (p1 - p2).length();
    const auto len_wi = (p3 - p2).length();
    const auto len_uw = (p1 - p3).length();

    auto max_value = std::acos(std::clamp<double>(
        dot((p3 - p2), (p1 - p2)) / len_ui /
        len_wi, -1, 1));
    auto radian = std::acos(std::clamp<double>(
        dot((p2 - p1), (p3 - p1)) / len_ui /
        len_uw, -1, 1));
    if (radian > max_value)
        max_value = radian;
    radian = std::acos(std::clamp<double>(
        dot((p1 - p3), (p2 - p3)) / len_uw /
        len_wi, -1, 1));
    if (radian > max_value)
        max_value = radian;

    if (max_value > 175. / 180. * M_PI ||
        mst.find_edge(triangle[0], triangle[2]).value().ref_time == 2 ||
        mst.find_edge(triangle[1], triangle[2]).value().ref_time == 2)
        return true;
    else
        return false;
}

void add_face(RSGraph& graph, const FaceType& item,
              std::vector<size_t>& faces)
{
    const auto& v_i = item[0];
    const auto& v_u = item[1];
    const auto& v_w = item[2];

    graph.increment_ref_time(v_i, v_w);
    graph.increment_ref_time(v_u, v_i);
    graph.increment_ref_time(v_w, v_u);
    faces.push_back(v_i);
    faces.push_back(v_u);
    faces.push_back(v_w);
}

struct RegisterFaceResult {
    InplaceVector<FaceType, 2> faces;
};
auto register_face(const RSGraph& mst, const NodeID v1, const NodeID v2) -> std::optional<RegisterFaceResult>
{
    GEL_ASSERT_EQ(mst.find_edge(v1, v2), RSGraph::InvalidEdgeID);
    auto shared_neighbors = mst.shared_neighbors_lazy(v1, v2);
    if (shared_neighbors.empty()) {
        return RegisterFaceResult(); // true
    }

    const Point& p1 = mst.m_vertices[v1].coords;
    const Point& p2 = mst.m_vertices[v2].coords;

    const auto possible_root1 = mst.predecessor(v1, v2).v;
    const auto angle1 = cal_radians_3d(p1 - mst.m_vertices[possible_root1].coords,
                                       mst.m_vertices[possible_root1].normal,
                                       p2 - mst.m_vertices[possible_root1].coords);
    const auto possible_root2 = mst.predecessor(v2, v1).v;
    const auto angle2 = cal_radians_3d(p2 - mst.m_vertices[possible_root2].coords,
                                       mst.m_vertices[possible_root2].normal,
                                       p1 - mst.m_vertices[possible_root2].coords);

    InplaceVector<FaceType, 2> temp;
    for (const auto v3 : shared_neighbors) {
        FaceType triangle{v1, v2, v3};
        if (v3 == possible_root1 && angle1 < M_PI) {
            if (routine_check(mst, triangle) ||
                mst.successor(v2, v1).v != v3 ||
                mst.successor(possible_root1, v2).v != v1) {
                return std::nullopt; // false
            }
            temp.emplace_back(FaceType{v1, v3, v2});
        }

        if (v3 == possible_root2 && angle2 < M_PI) {
            if (routine_check(mst, triangle) ||
                mst.successor(v1, v2).v != v3 ||
                mst.successor(possible_root2, v1).v != v2) {
                return std::nullopt; // false
            }
            temp.emplace_back(FaceType{v1, v2, v3});
        }
    }

    if (temp.empty())
        return std::nullopt; // false

    return RegisterFaceResult {
        temp,
    }; // true
}

/// @brief Connect handle to raise the genus number
/// @param smoothed_v: smoothed vertices of the point cloud
/// @param mst: graph and vertex information
/// @param kd_tree: kd-tree for knn query
/// @param connected_handle_root: [OUT] log the connected handles
/// @param k: number of kNN search
/// @param isEuclidean: if to use Euclidean distance
/// @param expected_genus expected genus of the mesh
/// @param step_thresh: step threshold for shortest distance path early stop
void connect_handle(
    const std::vector<Point>& smoothed_v,
    const Tree& kd_tree,
    RSGraph& mst,
    std::vector<NodeID>& connected_handle_root,
    const int k,
    const int step_thresh,
    const bool isEuclidean,
    const int expected_genus)
{
    std::vector<NodeID> imp_node;
    size_t num = 0;
    size_t edge_num = 0;
    // Collect vertices w/ an open angle larger than pi
    for (int i = 0; i < mst.no_nodes(); i++) {
        const auto& neighbors = mst.m_neighbors[i];
        // FIXME: Clean this up
        float last_angle = (*(--neighbors.end())).angle;
        float this_angle;

        for (auto& neighbor : neighbors) {
            this_angle = neighbor.angle;
            float angle_diff = this_angle - last_angle;
            if (angle_diff < 0)
                angle_diff += 2 * M_PI;
            if (angle_diff > M_PI)
                imp_node.push_back(i);
            last_angle = this_angle;
        }
        //for (const auto [neighbor_old, neighbor_next] : slide2_wraparound(neighbors)) {
        //    auto last_angle = neighbor_old.angle;
        //    auto this_angle = neighbor_next.angle;
        //    auto angle_diff = this_angle - last_angle;
        //    if (angle_diff < 0)
        //        angle_diff += 2 * M_PI;
        //    if (angle_diff > M_PI)
        //        imp_node.push_back(i);
        //}
    }

    struct Connection {
        NodeID p_from;
        NodeID p_to;
        uint tree_from;
        uint tree_to;
    };
    std::vector<Connection> connect;

    ThreadPool pool;
    const auto neighbors_map = calculate_neighbors(pool, smoothed_v, imp_node, kd_tree, k);

    GEL_ASSERT_EQ(imp_node.size(), neighbors_map.size());
    for (auto i = 0UL; i < imp_node.size(); i++) {
        const auto& this_v = imp_node[i];
        const auto& neighbors = neighbors_map[i];

        // Potential handle collection
        for (auto& neighbor: neighbors | std::views::drop(1)) {

            if (mst.find_edge(this_v, neighbor.id) != RSGraph::InvalidEdgeID)
                continue;
            const auto tree = mst.etf.representative(mst.predecessor(this_v, neighbor.id).tree_id);
            const auto to_tree = mst.etf.representative(mst.predecessor(neighbor.id, this_v).tree_id);
            const TEdge candidate(this_v, neighbor.id);
            if (tree != to_tree && geometry_check(mst, candidate, kd_tree)) {
                // TODO: Check if any tree shares root, and return corresponding edges
                connect.emplace_back(this_v, neighbor.id, tree, to_tree);
                break;
            }
        }
    }
    std::cout << "connections: " << connect.size() << "\n";

    // Select one handle
    Map<FaceConnectionKey, std::vector<int>, FaceConnectionKeyHasher> face_connections;
    for (int i = 0; i < connect.size(); i++) {
        uint tree = connect[i].tree_from; //tree_id[i];
        uint to_tree = connect[i].tree_to; //to_tree_id[i];
        if (to_tree > tree)
            std::swap(tree, to_tree);
        const auto key = FaceConnectionKey{tree, to_tree};
        if (!face_connections.contains(key))
            face_connections[key] = std::vector<int>{i};
        else {
            face_connections[key].push_back(i);
        }
    }

    // Sort
    std::vector<FacePair> sorted_face;
    for (const auto& key : face_connections | std::views::keys) {
        auto length = face_connections[key].size();
        sorted_face.emplace_back(length, key);
    }
    std::ranges::sort(sorted_face, face_comparator);
    for (const auto& val : sorted_face | std::views::values) {
        const std::vector<int>& idx_vec = face_connections[val];
        if (idx_vec.size() <= 5 || (expected_genus >= 0 && num >= expected_genus))
            break;

        bool isFind = false;
        for (const int idx : idx_vec) {
            const NodeID& this_v = connect[idx].p_from;
            const Point& query = mst.m_vertices[this_v].coords;
            const NodeID& connected_neighbor = connect[idx].p_to;
            // TODO: path and step_thresh are both dead
            const int steps = find_shortest_path(mst, this_v, connected_neighbor);
            if (steps > step_thresh) {
                isFind = true;
                const auto candidate = TEdge(this_v, connected_neighbor);
                if (geometry_check(mst, candidate, kd_tree)) {

                    mst.add_edge(this_v, connected_neighbor);
                    connected_handle_root.push_back(this_v);
                    connected_handle_root.push_back(connected_neighbor);

                    edge_num++;
                }
            }
        }
        if (isFind) {
            num++;
        }
    }

    std::cout << "Handle Connection done :)" << std::endl;
    std::cout << std::to_string(num) << " pairs of faces are connected." << std::endl;
    std::cout << std::to_string(edge_num) << " edges are connected." << std::endl;
}

bool explore(
    const RSGraph& G,
    const NodeID i,
    FacePriorityQueue& queue,
    const std::vector<ConnectionLength>& length_thresh,
    const bool is_euclidean)
{
    const NodeID v_i = i;
    bool isFound = false;
    for (const auto& neighbor : G.m_neighbors[i]) {
        const NodeID v_u = neighbor.v;
        const NodeID v_w = G.successor(i, v_u).v;

        const Point& w_pos = G.m_vertices[v_w].coords;
        const Point& u_pos = G.m_vertices[v_u].coords;
        const Point& i_pos = G.m_vertices[v_i].coords;
        const Vec3& i_normal = G.m_vertices[v_i].normal;
        const Vec3& u_normal = G.m_vertices[v_u].normal;
        const Vec3& w_normal = G.m_vertices[v_w].normal;
        const auto angle = cal_radians_3d(w_pos - i_pos, i_normal,
                                          u_pos - i_pos);
        const bool isLargerThanPi = angle < M_PI;
        const FaceType face_vector{v_i, v_u, v_w};
        if (v_u != v_w &&
            isLargerThanPi &&
            G.find_edge(v_u, v_w) == RSGraph::InvalidEdgeID) {
            const auto score = (is_euclidean)
                                   ? (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length()
                                   : cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                                                   u_normal, w_normal);
            if (score > length_thresh[v_u].max_length || score > length_thresh[v_w].max_length)
                continue;
            else if (score >= 0) {
                queue.emplace(face_vector, score);
                isFound = true;
            }
        }
    }

    return isFound;
}

bool check_branch_validity(const RSGraph& G, const NodeID root, const NodeID branch1, const NodeID branch2)
{
    constexpr auto angle_thresh = 0. / 180. * M_PI;

    const Point& pos_u = G.m_vertices[branch1].coords;
    const Point& pos_w = G.m_vertices[branch2].coords;
    const Vec3& normal_u = G.m_vertices[branch1].normal;
    const Vec3& normal_w = G.m_vertices[branch2].normal;

    const auto is_valid = [&root](auto this_radian, auto former, auto next) {
        if (next.v == root) {
            auto diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                return true;
        }
        if (former.v == root) {
            auto diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                return true;
        }
        return false;
    };

    // Check u's RS validity
    {
        const auto this_radian = cal_radians_3d(pos_w - pos_u, normal_u);
        const auto& former = G.predecessor(branch1, branch2);
        const auto& next = G.successor(branch1, branch2);
        if (!is_valid(this_radian, former, next))
            return false;

        // Thresh on angle
        auto diff_angle_thresh = this_radian - former.angle;
        if (diff_angle_thresh < 0)
            diff_angle_thresh += M_PI * 2.;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }

    // Check w's RS validity
    {
        const auto this_radian = cal_radians_3d(pos_u - pos_w, normal_w);
        const auto& former = G.predecessor(branch2, branch1);
        const auto& next = G.successor(branch2, branch1);
        if (!is_valid(this_radian, former, next))
            return false;

        // Thresh on angle
        if (const auto diff_angle_thresh2 = -this_radian + next.angle + 2. * M_PI;
            diff_angle_thresh2 < angle_thresh)
            return false;
    }

    return true;
}


// TODO: validity of what?
bool check_validity(
    const RSGraph& graph,
    const std::pair<FaceType, float>& item)
{
    const NodeID v_i = item.first[0];
    const NodeID v_u = item.first[1];
    const NodeID v_w = item.first[2];
    const Point& pos_i = graph.m_vertices[v_i].coords;
    const Point& pos_u = graph.m_vertices[v_u].coords;
    const Point& pos_w = graph.m_vertices[v_w].coords;
    const Vec3& normal_i = graph.m_vertices[v_i].normal;

    if (graph.find_edge(v_u, v_w) != RSGraph::InvalidEdgeID)
        return false;

    // Non-manifold edge check
    if (graph.find_edge(v_i, v_u).value().ref_time == 2 ||
        graph.find_edge(v_i, v_w).value().ref_time == 2)
        return false;

    // Check this rotation system
    if (graph.successor(v_i, v_u).v != v_w)
        return false;
    const auto angle = cal_radians_3d(pos_w - pos_i, normal_i, pos_u - pos_i);
    if (angle > M_PI)
        return false;

    // Check the rotation system's validity of branch nodes
    if (!check_branch_validity(graph, v_i, v_u, v_w)) {
        return false;
    }

    return true;
}

void triangulate(
    std::vector<size_t>& faces,
    RSGraph& graph,
    const bool is_euclidean,
    const std::vector<ConnectionLength>& length_thresh,
    const std::vector<NodeID>& connected_handle_root)
{
    // Since we are clearing every element here, a TreeSet has much better performance than a Swiss table
    OrderedSet<NodeID> to_visit;
    FacePriorityQueue queue;

    // Init priority queue
    for (auto i = 0UL; i < graph.no_nodes(); i++) {
        to_visit.insert(i);
    }

    for (const auto i : connected_handle_root) {
        explore(graph, i, queue, length_thresh, is_euclidean);
        to_visit.erase(i);
    }

    std::cout << "Global init done :)" << std::endl;

    while (!to_visit.empty()) {
        while (!queue.empty()) {
            std::pair<FaceType, float> item = queue.top();
            queue.pop();

            if (item.second >= 0.0 && !check_validity(graph, item)) {
                continue;
            }

            // Add the edge
            const NodeID v_i = item.first[0];
            const NodeID v_u = item.first[1];
            const NodeID v_w = item.first[2];

            if (graph.find_edge(v_u, v_w) == RSGraph::InvalidEdgeID) {
                graph.add_edge(v_u, v_w);
                add_face(graph, item.first, faces);
            } else {
                continue;
            }

            // Deal with incident triangles
            for (const NodeID incident_root : graph.shared_neighbors_lazy(v_u, v_w)) {
                if (incident_root == v_i)
                    continue;
                std::array face{incident_root, v_w, v_u};

                // Non-manifold edge check
                if (graph.find_edge(incident_root, v_u).value().ref_time == 2 ||
                    graph.find_edge(incident_root, v_w).value().ref_time == 2 ||
                    graph.find_edge(v_u, v_w).value().ref_time == 2)
                    continue;

                add_face(graph, face, faces);
            }

            to_visit.erase(v_u);
            to_visit.erase(v_w);

            // Explore and sanity check
            explore(graph, v_u, queue, length_thresh, is_euclidean);
            explore(graph, v_w, queue, length_thresh, is_euclidean);
        }

        if (!to_visit.empty()) {
            NodeID pick = *to_visit.begin();
            to_visit.erase(pick);
            explore(graph, pick, queue, length_thresh, is_euclidean);
        }
    }
}

/// @brief Build minimum spanning tree (MST)
/// @param out_mst: [OUT] constructed MST
/// @param g: connection information of the mst
/// @param root root node
/// @param is_euclidean: if to use Euclidean distance
/// @param vertices: coordinates of the point cloud
/// @param normals: normal of the point cloud
void build_mst(
    const SimpGraph& g,
    const NodeID root,
    RSGraph& out_mst,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    const bool is_euclidean)
{
    RSGraph temp = minimum_spanning_tree(g, root, normals, vertices, is_euclidean);

    // Fix strong ambiguous points
    if (!is_euclidean) {
        for (int i = 0; i < temp.m_edges.size(); i++) {
            const NodeID& source = temp.m_edges[i].source;
            const NodeID& target = temp.m_edges[i].target;
            const Vec3& normal1 = temp.m_vertices[source].normal;
            const Vec3& normal2 = temp.m_vertices[target].normal;
            const Point& pos1 = temp.m_vertices[source].coords;
            const Point& pos2 = temp.m_vertices[target].coords;
            if (temp.valence(source) >= 2 && temp.valence(target) >= 2)
                continue;
            const Vec3 edge = pos2 - pos1;

            const Vec3 normal_sum = (normal1 + normal2) * 0.5;
            const auto cos_angle = std::abs(dot(edge, normal_sum / edge.length()));
            if (cos_angle > std::cos(10. / 180. * M_PI)) {
                NodeID parent;
                if (temp.valence(source) == 1) {
                    temp.m_vertices[source].normal = temp.m_vertices[target].normal;
                    parent = target;
                } else {
                    temp.m_vertices[target].normal = temp.m_vertices[source].normal;
                    parent = source;
                }

                for (const auto neighbor : temp.neighbors_lazy(parent)) {
                    if (temp.m_vertices[neighbor].normal_rep == Vertex::InvalidNormalRep) {
                        // this conversion is fine since we need 10^18 vertices for it to matter
                        temp.m_vertices[neighbor].normal_rep = parent;
                    } else {
                        // Collision!
                        temp.m_vertices[neighbor].normal_rep = Vertex::CollisionRep;
                    }
                }
            }
        }
        for (int i = 0; i < temp.m_vertices.size(); i++) {
            if (temp.m_vertices[i].normal_rep >= 0)
                temp.m_vertices[i].normal = temp.m_vertices[temp.m_vertices[i].normal_rep].normal;
        }
    }

    // Build corrected MST
    size_t verts = 0;
    for (const auto& vertex : temp.m_vertices) {
        out_mst.add_node(vertex.coords, vertex.normal);
        verts++;
    }

    size_t edges = 0;
    for (const auto& m_edge : temp.m_edges) {
        out_mst.add_edge(m_edge.source, m_edge.target);
        edges++;
    }
    GEL_ASSERT_EQ(verts, edges + 1);
}


auto estimate_normals_included_normals(
    const std::vector<Point>& vertices,
    std::vector<Vec3>& normals,
    const Distance dist) -> std::vector<Point>
{
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    ThreadPool pool;
    std::vector<Vec3> smoothed_v;
    std::vector<Point> temp;
    const int smoothing_size = std::max(static_cast<int>(static_cast<double>(vertices.size()) / 2000.), 192);

    normalize_normals(normals);
    if (dist == Distance::Euclidean) {
        smoothed_v = vertices;
    } else if (dist == Distance::Tangent) {
        const auto indices = std::ranges::iota_view(0UL, vertices.size());
        const auto kdTree = build_kd_tree_of_indices(vertices, indices);
        const auto neighbors =
            calculate_neighbors(pool, vertices, kdTree, smoothing_size);
        weighted_smooth(pool, vertices, normals, neighbors, smoothed_v);


        temp.reserve(smoothed_v.size());
        std::swap(temp, smoothed_v);
        smoothed_v.clear();

        weighted_smooth(pool, temp, normals, neighbors, smoothed_v);
    } else {
        GEL_ASSERT(false, "unreachable");
    }
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    return smoothed_v;
}

auto estimate_normals_no_normals(
    const std::vector<Point>& vertices,
    std::vector<Vec3>& normals,
    const Distance dist) -> std::vector<Point>
{
    GEL_ASSERT_EQ(normals.size(), 0LU);
    ThreadPool pool;
    const int smoothing_size = std::max(static_cast<int>(static_cast<double>(vertices.size()) / 2000.), 192);

    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    const auto kdTree = build_kd_tree_of_indices(vertices, indices);
    auto neighbors = calculate_neighbors(pool, vertices, kdTree, smoothing_size);

    estimate_normal_no_normals_memoized(pool, vertices, neighbors, normals);

    std::vector<Point> smoothed_v;

    switch (dist) {
    case Distance::Euclidean:
        smoothed_v = vertices;
        break;
    case Distance::Tangent:
        std::cout << "Smoothing round 1 ..." << std::endl;
        weighted_smooth(pool, vertices, normals, neighbors, smoothed_v);
        break;
    default:
        GEL_ASSERT(false, "unreachable");
    }

    const auto temp_tree1 = build_kd_tree_of_indices(smoothed_v, indices);
    neighbors = calculate_neighbors(pool, smoothed_v, temp_tree1, smoothing_size);
    estimate_normal_no_normals_memoized(pool, smoothed_v, neighbors, normals);

    if (dist == Distance::Tangent) {
        std::vector<Point> temp;
        temp.reserve(smoothed_v.size());
        std::swap(temp, smoothed_v);
        smoothed_v.clear();
        std::cout << "Smoothing round 2 ..." << std::endl;
        weighted_smooth(pool, temp, normals, neighbors, smoothed_v);

        const Tree temp_tree2 = build_kd_tree_of_indices(smoothed_v, indices);
        neighbors = calculate_neighbors(pool, smoothed_v, temp_tree2, smoothing_size, std::move(neighbors));
        estimate_normal_no_normals_memoized(pool, smoothed_v, neighbors, normals);
    }
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    return smoothed_v;
}

[[nodiscard]]
auto estimate_normals_and_smooth(
    IExecutor& pool,
    const std::vector<Point>& org_vertices,
    std::vector<Vec3>& org_normals,
    Distance dist) -> std::vector<Point>
{
    if (org_normals.empty()) {
        return estimate_normals_no_normals(org_vertices, org_normals, dist);
    } else {
        return estimate_normals_included_normals(org_vertices, org_normals, dist);
    }
}

struct Components {
    std::vector<std::vector<Point>> vertices;
    std::vector<std::vector<Point>> smoothed_v;
    std::vector<std::vector<Vec3>> normals;
};

/// @brief Find the number of connected components and separate them
/// @param pool: Thread pool to use
/// @param vertices: vertices of the point cloud
/// @param smoothed_v: smoothed vertices of the point cloud
/// @param normals: normal of the point cloud vertices
/// @param neighbor_map
/// @param opts: theta: (cross-connection threshold) angle threshold to avoid connecting vertices on different surface
///              r: (outlier_thresh) threshold distance(?) to remove an outlier
///              k
///              isEuclidean
/// @return Split components
[[nodiscard]]
auto split_components(
    IExecutor& pool,
    const NeighborMap& neighbor_map,
    std::vector<Point>&& vertices,
    std::vector<Vec3>&& normals,
    std::vector<Point>&& smoothed_v,
    const RsROpts& opts)
    -> Components
{
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    GEL_ASSERT_EQ(vertices.size(), smoothed_v.size());
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<std::vector<Vec3>> component_normals;

    // Identifies clusters of vertices which are reconstructed to disparate meshes
    const double outlier_thresh = opts.r;
    double total_edge_length = 0;
    AMGraph::NodeSet sets;
    SimpGraph components;
    for (int i = 0; i < vertices.size(); i++) {
        sets.insert(components.add_node());
    }

    // Construct graph
    for (const auto& neighbors : neighbor_map) {
        const NodeID this_idx = neighbors[0].id;
        for (const auto& neighbor : neighbors | std::views::drop(1)) {
            const NodeID idx = neighbor.id;
            if (this_idx < idx) continue;
            const double length = neighbor.distance;

            total_edge_length += length;

            components.connect_nodes(this_idx, idx);
        }
    }
    const double thresh_r = (total_edge_length * outlier_thresh) / static_cast<double>(components.no_edges());
    // Remove Edges Longer than the threshold
    std::vector<std::pair<NodeID, NodeID>> edge_rm_v;
    for (NodeID vertex1 = 0; vertex1 < components.no_nodes(); ++vertex1) {
        for (const auto& vertex2 : components.neighbors_lazy(vertex1)) {
            const double edge_length = (vertices[vertex1] - vertices[vertex2]).length();

            if (edge_length > thresh_r) {
                edge_rm_v.emplace_back(vertex1, vertex2);
            }
        }
    }
    for (auto& [fst, snd] : edge_rm_v) {
        components.disconnect_nodes(fst, snd);
    }
    // Find Components
    const auto components_vec = connected_components(components.inner(), sets);
    std::cout << "The input contains " << components_vec.size() << " connected components." << std::endl;
    // Valid Components and create new vectors for components
    const auto threshold = std::min<size_t>(vertices.size(), 100);
    for (auto& component : components_vec) {
        if (component.size() >= threshold) {
            std::vector<Vec3> this_normals;
            this_normals.reserve(component.size());
            std::vector<Point> this_vertices;
            this_vertices.reserve(component.size());
            std::vector<Point> this_smoothed_v;
            this_smoothed_v.reserve(component.size());
            for (const auto& element : component) {
                this_vertices.push_back(vertices[element]);
                this_smoothed_v.push_back(smoothed_v[element]);
                this_normals.push_back(normals[element]);
            }
            component_normals.emplace_back(std::move(this_normals));
            component_vertices.emplace_back(std::move(this_vertices));
            component_smoothed_v.emplace_back(std::move(this_smoothed_v));
        }
    }
    std::cout << component_vertices.size() << " of them will be reconstructed." << std::endl;

    return Components(
        std::move(component_vertices),
        std::move(component_smoothed_v),
        std::move(component_normals)
    );
}

void export_edges(const RSGraph& g, const std::string& out_path) {
    std::ofstream file(out_path);
    // Write vertices
    file << "# List of geometric vertices" << std::endl;
    for (int i = 0; i < g.m_vertices.size(); i++) {
        Point this_coords = g.m_vertices[i].coords;
        file << "v " << std::to_string(this_coords[0])
            << " " << std::to_string(this_coords[1])
            << " " << std::to_string(this_coords[2]) << "\n";
    }

    // Write lines
    file << std::endl;
    file << "# Line element" << std::endl;
    for (int i = 0; i < g.m_vertices.size(); i += 2)
        file << "l " << i + 1
        << " " << i + 2 << "\n";
    //file.close();
}

void export_graph(const RSGraph& g, const std::string& out_path) {
    std::ofstream file(out_path);
    // Write vertices
    file << "# List of geometric vertices\n";
    for (int i = 0; i < g.m_vertices.size(); i++) {
        Point this_coords = g.m_vertices[i].coords;
        file << "v " << std::to_string(this_coords[0])
            << " " << std::to_string(this_coords[1])
            << " " << std::to_string(this_coords[2]) << "\n";
    }

    // Write lines
    file << "\n# Line elements\n";
    for (int i = 0; i < g.m_edges.size(); ++i) {
        auto edge = g.m_edges[i];
        file << "l " << (edge.source + 1) << " " << (edge.target + 1) << "\n";
    }
}

bool vec3_eq_(const Vec3& lhs, const Vec3& rhs, double eps = 1e-4)
{
    return lhs.all_le(rhs + Vec3(eps)) && lhs.all_ge(rhs - Vec3(eps));
}


RSGraph from_simp_graph(const SimpGraph& graph, const std::vector<Point>& points)
{
    RSGraph copy;

    for (auto id: graph.node_ids()) {
        copy.add_node(points.at(id));
    }
    for (auto id: graph.node_ids()) {
        for (auto neighbor: graph.neighbors_lazy(id)) {
            if (id < neighbor) {
                copy.add_edge(id, neighbor);
            }
        }
    }
    return copy;
}

auto component_to_manifold(
    IExecutor& pool,
    const RsROpts& opts,
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& smoothed_v,
    const Tree& kd_tree,
    const NeighborMap& neighbor_map
) -> Manifold
{
    Util::RSRTimer inner_timer;
    std::cout << "Init mst" << std::endl;
    // Initial Structure
    RSGraph mst;
    mst.reserve(vertices.size(), opts.k);

    std::vector<PEdgeLength> edge_length;
    std::vector<ConnectionLength> connection_lengths(vertices.size(), ConnectionLength());
    {
        inner_timer.start("init_graph");
        SimpGraph g = init_graph(smoothed_v, normals, neighbor_map, connection_lengths, opts.k, opts.theta,
                                 opts.dist == Distance::Euclidean);
        inner_timer.end("init_graph");
        // Generate MST
        inner_timer.start("build_mst");

        build_mst(g, 0, mst, normals, smoothed_v, opts.dist == Distance::Euclidean);
        //export_graph(from_simp_graph(g, vertices), "all_connected.obj");
        //export_graph(mst, "mst.obj");
        inner_timer.end("build_mst");

        // Edge arrays and sort
        inner_timer.start("edge_length");
        for (NodeID node : g.node_ids()) {
            for (NodeID node_neighbor : g.neighbors_lazy(node)) {
                if (node < node_neighbor) {
                    const Vec3 edge = smoothed_v[node] - smoothed_v[node_neighbor];
                    const double len = (opts.dist == Distance::Euclidean)
                                           ? edge.length()
                                           : cal_proj_dist(edge, normals[node], normals[node_neighbor]);

                    if (len > connection_lengths[node].pre_max_length ||
                        len > connection_lengths[node_neighbor].pre_max_length)
                        continue;
                    edge_length.emplace_back(TEdge(node, node_neighbor), len);
                }
            }
        }
        inner_timer.end("edge_length");

        inner_timer.start("sort");
        std::ranges::sort(edge_length, edge_comparator);
        inner_timer.end("sort");
    }

    // Initialize face loop label
    mst.etf.reserve(6 * vertices.size());
    init_face_loop_label(mst);

    std::vector<size_t> flattened_face;
    flattened_face.reserve(vertices.size() * 2 * 3); // A fairly reasonable heuristic

    std::vector<std::pair<Point, NodeID>> neighbors;
    inner_timer.start("edge_connection");
     for (auto& [this_edge, edge_len] : edge_length) {
         if (mst.find_edge(this_edge.from, this_edge.to) == RSGraph::InvalidEdgeID &&
             vanilla_check(mst, this_edge, kd_tree, neighbors)) {
             auto result = register_face(mst, this_edge.from, this_edge.to);
             if (result.has_value()) {
                 mst.add_edge(this_edge.from, this_edge.to);
                 mst.maintain_face_loop(this_edge.from, this_edge.to);
                 for (const auto& face : result->faces) {
                     add_face(mst, face, flattened_face);
                 }
             }
         }
     }
    inner_timer.end("edge_connection");
    std::cout << "edge length length: " << edge_length.size() << "\n";

    // Create handles & Triangulation
    inner_timer.start("triangulation");
    if (opts.genus != 0) {
        std::vector<NodeID> connected_handle_root;
        connect_handle(smoothed_v, kd_tree, mst, connected_handle_root, opts.k, opts.n,opts.dist == Distance::Euclidean, opts.genus);
        triangulate(flattened_face, mst, opts.dist == Distance::Euclidean, connection_lengths, connected_handle_root);
    }
    inner_timer.end("triangulation");

    inner_timer.start("build_manifold");
    Manifold res;

    build_manifold(res, vertices, flattened_face, 3);
    inner_timer.end("build_manifold");

    std::cout << "\n";
    inner_timer.show();
    std::cout << "\n";

    return res;
}

auto point_cloud_to_mesh_impl(
    std::vector<Point>&& vertices_copy,
    std::vector<Vec3>&& normals_copy,
    Util::RSRTimer& timer,
    ThreadPool& pool,
    const RsROpts& opts) -> Manifold
{
    Manifold output;

    auto in_smoothed_v = estimate_normals_and_smooth(pool, vertices_copy, normals_copy, opts.dist);
    const Tree kd_tree = build_kd_tree_of_indices(in_smoothed_v, std::views::iota(0UL, in_smoothed_v.size()));

    // Find components
    timer.start("Split components");
    std::cout << "Split components\n";
    // Note: the cross connection filtering needs to be synced with the inner loop, else there are issues
    auto neighbor_map = calculate_neighbors(pool, in_smoothed_v, kd_tree, opts.k);
    Parallel::foreach(pool, std::views::iota(0UL, vertices_copy.size()), [&](const NodeID idx) {
        filter_out_cross_connection(neighbor_map[idx], normals_copy, idx, opts.theta, opts.dist == Distance::Euclidean);
    });

    auto [component_vertices,
          component_smoothed_v,
          component_normals] =
        split_components(pool,
                         neighbor_map,
                         std::move(vertices_copy),
                         std::move(normals_copy),
                         std::move(in_smoothed_v),
                         opts);
    timer.end("Split components");
    // There is no guarantee that there is more than one component, and components can
    // be highly non-uniform in terms of how many primitives they have. That means we cannot
    // rely on this loop for good parallelization opportunities.
    timer.start("Algorithm");
    for (size_t component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " << std::to_string(component_id) << " ... (" << component_vertices[component_id].size() << " vertices)" << std::endl;

        std::vector<Point> vertices_of_this = std::move(component_vertices[component_id]);
        std::vector<Vec3> normals_of_this = std::move(component_normals[component_id]);
        std::vector<Point> smoothed_v_of_this = std::move(component_smoothed_v[component_id]);
        GEL_ASSERT(vertices_of_this.size() == normals_of_this.size());
        GEL_ASSERT(vertices_of_this.size() == smoothed_v_of_this.size());

        // While I would like to move this up, there are some nontrivial changes made to each component
        // inside split_components eve if we only have one component to worry about.
        const auto indices_of_this = std::ranges::iota_view(0UL, smoothed_v_of_this.size());
        Tree kd_tree_of_this = build_kd_tree_of_indices(smoothed_v_of_this, indices_of_this);

        auto neighbor_map_of_this = calculate_neighbors(pool, smoothed_v_of_this, kd_tree_of_this, opts.k);
        //Parallel::foreach(pool, std::views::iota(0UL, smoothed_v_of_this.size()), [&](NodeID idx) {
        //    filter_out_cross_connection(neighbor_map_of_this[idx], normals_of_this, idx, opts.theta,
        //                                opts.dist == Distance::EUCLIDEAN);
        //});

        auto res = component_to_manifold(
            pool,
            opts,
            vertices_of_this,
            normals_of_this,
            smoothed_v_of_this,
            kd_tree_of_this,
            neighbor_map_of_this);

        output.merge(res);
        // if (component_id == 0) {
        //     std::swap(output, res);
        // } else {
        //     output.merge(res);
        // }
    }
    timer.end("Algorithm");

    return output;
}

auto point_cloud_to_mesh(
    const std::vector<Point>& vertices_in,
    const std::vector<Vec3>& normals_in,
    const RsROpts& opts) -> Manifold
{
    Util::RSRTimer timer;
    ThreadPool pool;
    timer.start("Whole process");
    timer.start("Validation");
    if (!normals_in.empty()) {
        GEL_ASSERT_EQ(vertices_in.size(), normals_in.size(), "Vertices and normals must be the same size");
    }
    for (const auto& point: vertices_in) {
        GEL_ASSERT_FALSE(std::isnan(point[0]) || std::isnan(point[1]) || std::isnan(point[2]), "Bad point input");
    }
    for (const auto& normal: normals_in) {
        GEL_ASSERT_FALSE(std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2]), "Bad normal input");
        GEL_ASSERT_NEQ(normal, Vec3(0.0), "Bad normal input");
    }
    auto vertices_copy = vertices_in;
    auto normals_copy = normals_in;
    timer.end("Validation");


    // Estimate normals & orientation & weighted smoothing
    timer.start("Estimate and smooth normals");
    std::vector<Point> in_smoothed_v = estimate_normals_and_smooth(pool, vertices_copy, normals_copy, opts.dist);
    timer.end("Estimate and smooth normals");


    const auto indices = std::ranges::iota_view(0UL, in_smoothed_v.size());
    const Tree kd_tree = build_kd_tree_of_indices(in_smoothed_v, indices);
    if (normals_in.empty()) {
        std::cout << "correct normal orientation\n";
        timer.start("Correct normal orientation");
        correct_normal_orientation(pool, kd_tree, in_smoothed_v, normals_copy, opts.k);
        timer.end("Correct normal orientation");
    }

    Manifold result = point_cloud_to_mesh_impl(
        std::move(vertices_copy),
        std::move(normals_copy),
        timer,
        pool,
        opts);
    timer.end("Whole process");
    const std::string line(40, '=');
    std::cout << line << "\n\n";
    timer.show();
    return result;
}

template <typename Collection, typename IndexRange>
auto indexed_select(const Collection& collection,
                    const IndexRange& indices) -> std::vector<typename Collection::value_type>
{
    std::vector<typename Collection::value_type> result;
    result.reserve(collection.size());
    for (auto idx : indices) {
        result.push_back(collection.at(idx));
    }
    return result;
}

std::vector<Vec3> validate_normals(ThreadPool& pool, const std::vector<Point>& vertices, const std::vector<Vec3>& normals, const RsROpts& reconstruction_options)
{
    auto check = [&]() {
        if (!normals.empty() && vertices.size() != normals.size()) {
            std::cerr << "Bad normal input: normals must either be empty or the same size as vertices" << "\n";
            return false;
        }
        for (const auto& point: vertices) {
            GEL_ASSERT_FALSE(std::isnan(point[0]) || std::isnan(point[1]) || std::isnan(point[2]), "Bad point input: Found a NaN coordinate.");
        }
        for (const auto& normal: normals) {
            if (std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2])) {
                std::cerr << "Bad normal input: normal is NaN" << "\n";
                return false;
            }
            if (normal == Vec3(0.0)) {
                std::cerr << "Bad normal input: normal is 0 vector" << "\n";
                return false;
            }
        }
        return true;
    };
    if (check()) {
        return normals;
    } else {
        std::cout << "generating normals...";
        auto normals_copy = std::vector<Vec3>();
        auto smoothed_points = estimate_normals_and_smooth(pool, vertices, normals_copy, reconstruction_options.dist);
        return normals_copy;
    }

}

auto point_cloud_collapse_reexpand(
    const std::vector<Point>& vertices_in,
    const std::vector<Vec3>& normals_in,
    const CollapseOpts& collapse_options,
    const RsROpts& rsr_opts,
    const ReexpandOptions& reexpand_opts) -> Manifold
{
    if (collapse_options.max_iterations == 0) {
        return point_cloud_to_mesh(vertices_in, normals_in, rsr_opts);
    }
    Util::RSRTimer timer;
    ThreadPool pool;

    timer.start("Whole process");
    timer.start("Validation");
    if (!normals_in.empty()) {
        GEL_ASSERT_EQ(vertices_in.size(), normals_in.size(), "Vertices and normals must be the same size");
    }
    for (const auto& point: vertices_in) {
        GEL_ASSERT_FALSE(std::isnan(point[0]) || std::isnan(point[1]) || std::isnan(point[2]), "Bad point input");
    }
    for (const auto& normal: normals_in) {
        GEL_ASSERT_FALSE(std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2]), "Bad normal input");
        GEL_ASSERT_NEQ(normal, Vec3(0.0), "Bad normal input");
    }
    auto vertices_copy = vertices_in;
    auto normals_copy = normals_in;
    timer.end("Validation");


    // Estimate normals & orientation & weighted smoothing
    timer.start("Estimate and smooth normals");
    std::vector<Point> in_smoothed_v = estimate_normals_and_smooth(pool, vertices_copy, normals_copy, collapse_options.distance);
    timer.end("Estimate and smooth normals");



    if (normals_in.empty()) {
        std::cout << "correct normal orientation\n";
        timer.start("Correct normal orientation");
        const auto indices = std::ranges::iota_view(0UL, in_smoothed_v.size());
        const Tree kd_tree = build_kd_tree_of_indices(in_smoothed_v, indices);
        correct_normal_orientation(pool, kd_tree, in_smoothed_v, normals_copy, rsr_opts.k);
        timer.end("Correct normal orientation");
    }

    timer.start("Collapse");
    auto [collapse, point_cloud] = collapse_points(vertices_copy, normals_copy, collapse_options);
    timer.end("Collapse");

    auto [points_collapsed, normals_collapsed, indices_collapsed] = std::move(point_cloud);

    Manifold manifold = point_cloud_to_mesh_impl(
        std::move(points_collapsed),
        std::move(normals_collapsed),
        timer,
        pool,
        rsr_opts);

    timer.start("Reexpand");
    if (reexpand_opts.enabled)
        reexpand_points(manifold, std::move(collapse), reexpand_opts);
    timer.end("Reexpand");

    timer.end("Whole process");
    const std::string line(40, '=');
    std::cout << line << "\n\n";
    timer.show();
    return manifold;
}

auto point_cloud_normal_estimate(const std::vector<Point>& vertices,
                                 const std::vector<Vec3>& normals,
                                 const bool is_euclidean) -> NormalEstimationResult
{
    auto normals_copy = normals;
    const auto dist = is_euclidean ? Distance::Euclidean : Distance::Tangent;
    ThreadPool pool;
    const auto smoothed_v = estimate_normals_and_smooth(pool, vertices, normals_copy, dist);
    return {vertices, normals_copy, smoothed_v};
}
} // namespace GEL::HMesh::RsR
