#include <GEL/Util/RSRTimer.h>
#include <GEL/Geometry/normal.h>
#include <GEL/HMesh/RsR2.h>
#include <GEL/Util/ParallelAdapters.h>
#include <ranges> // std::views

namespace HMesh::RSR
{
// Face Loop
const Neighbor& successor(const RSGraph& g, const NodeID& root, const NodeID& branch);
const Neighbor& predecessor(const RSGraph& g, const NodeID& root, const NodeID& branch);
void maintain_face_loop(RSGraph& g, NodeID source, NodeID target);
const Neighbor& get_neighbor_info(const RSGraph& g, const NodeID& root, const NodeID& branch);
void find_common_neighbor(const RSGraph& g, NodeID neighbor, NodeID root, std::vector<NodeID>& shared_neighbors);

struct FaceComparator {
    bool operator()(const std::pair<FaceType, float>& left,
                    const std::pair<FaceType, float>& right) const
    {
        return (left.second) > (right.second);
    }
};

using FacePriorityQueue = std::priority_queue<
    std::pair<FaceType, float>,
    std::vector<std::pair<FaceType, float>>,
    FaceComparator>;

using PEdgeLength = std::pair<float, int>;

struct FaceConnectionKey {
    uint tree_id;
    uint to_tree_id;

    bool operator==(const FaceConnectionKey& other) const noexcept
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

using namespace ::HMesh;
using ::Util::AttribVec;
using ::HMesh::Manifold;

inline bool edge_comparator(const PEdgeLength& l, const PEdgeLength& r)
{
    return l.first < r.first;
}

inline bool face_comparator(const FacePair& l, const FacePair& r)
{
    return l.first > r.first;
}

inline bool neighbor_comparator(const NeighborPair& l, const NeighborPair& r)
{
    return l.first > r.first;
}

/**
    * @brief Calculate the reference vector for the rotation system
    *
    * @param normal: normal direction for the target vertex
    *
    * @return the reference vector
    */
Vec3 calculate_ref_vec(const Vec3& normal)
{
    constexpr float eps = 1e-6;
    double second = normal[1];
    if (second == 0.)
        second += eps;
    auto ref_vec = Vec3(0, -normal[2] / second, 1);
    if (normal[2] == 1.)
        ref_vec = Vec3(0., 1., 0.);
    ref_vec.normalize();
    return ref_vec;
}

/**
    * @brief Calculate the radian in the rotation system
    *
    * @param branch: vector of the out-going edge
    * @param normal: normal of the root vertex
    *
    * @return radian
    */
double cal_radians_3d(const Vec3& branch, const Vec3& normal)
{
    const Vec3 proj_vec = branch - dot(normal, branch) /
        normal.length() * normal;

    const auto ref_vec = calculate_ref_vec(normal);

    if (proj_vec.length() == 0.0)
        return 0.;

    const Vec3 proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    const auto value = std::clamp<double>(dot(proj_vec, proj_ref) / proj_vec.length() /
                                          proj_ref.length(), -1, 1);
    double radian = std::acos(value);
    if (dot(cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;

    if (std::isnan(radian)) [[unlikely]] {
        std::cout << normal << std::endl;
        std::cout << ref_vec << std::endl;
        std::cout << "error" << std::endl;
    }
    return radian;
}

/**
    * @brief Calculate the radian given the reference vector
    *
    * @param branch_vec: vector of the out-going edge
    * @param normal: normal of the root vertex
    * @param ref_vec: the reference vector
    *
    * @return radian
    */
double cal_radians_3d(const Vec3& branch_vec, const Vec3& normal, const Vec3& ref_vec)
{
    const Vec3 proj_vec = branch_vec - dot(normal, branch_vec) /
        normal.length() * normal;
    if (std::abs(proj_vec.length()) < 1e-8)
        return 0.;

    const Vec3 proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    const auto value = std::clamp<double>(
        dot(proj_vec, proj_ref) / proj_vec.length() /
        proj_ref.length(), -1, 1);
    double radian = std::acos(value);
    if (dot(CGLA::cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;
    return radian;
}


/**
    * @brief neighbor search within a specific radius
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param dist: the radius of the search ball
    * @param neighbors: [OUT] indices of k nearest neighbors
    *
    * @return None
    */
void nn_search(const Point& query, const Tree& kdTree,
               const double dist, std::vector<NodeID>& neighbors)
{
    std::vector<Point> neighbor_coords;
    std::vector<double> neighbor_dist;
    kdTree.in_sphere(query, dist, neighbor_coords, neighbors);

    std::vector<NeighborPair> paired;
    auto i = 0UL;
    for (auto& neighbor_coord : neighbor_coords) {
        neighbor_dist.push_back((neighbor_coord - query).length());
        paired.emplace_back(neighbor_dist[i], neighbors[i]);
        i++;
    }

    std::ranges::sort(paired, neighbor_comparator);
    for (size_t idx = 0; idx < paired.size(); ++idx) {
        neighbor_dist[idx] = paired[idx].first;
        neighbors[idx] = paired[idx].second;
    }
}

/**
    * @brief Calculate projection distance
    *
    * @param edge: the edge to be considered
    * @param this_normal: normal of one vertex
    * @param neighbor_normal: normal of another vertex
    *
    * @return projection distance
    */
double cal_proj_dist(const Vec3& edge, const Vec3& this_normal, const Vec3& neighbor_normal)
{
    const double Euclidean_dist = edge.length();
    const double neighbor_normal_length = dot(edge, normalize(neighbor_normal));
    const double normal_length = dot(edge, normalize(this_normal));
    double projection_dist = sqrt((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
    projection_dist += sqrt((Euclidean_dist * Euclidean_dist) -
        (neighbor_normal_length * neighbor_normal_length));
    projection_dist /= 2.;
    if (std::abs(dot(normalize(this_normal), normalize(neighbor_normal))) < std::cos(15. / 180. * M_PI))
        projection_dist = Euclidean_dist;
    return projection_dist;
}

auto filter_out_cross_connection(
    NeighborArray& neighbors,
    const std::vector<Vec3d>& normals,
    const NodeID this_idx,
    const double cross_conn_thresh,
    const bool isEuclidean
)
{
    const Vec3 this_normal = normals[this_idx];
    std::erase_if(neighbors, [&](auto& neighbor_info) {
        const Vec3 neighbor_normal = normals[neighbor_info.id];
        const double cos_theta = dot(this_normal, neighbor_normal) /
            this_normal.length() / neighbor_normal.length();
        const double cos_thresh =
            (isEuclidean) ? 0.0 : cos(cross_conn_thresh / 180. * M_PI);
        return cos_theta < cos_thresh;
    });
}


/**
    * @brief initialize the graph and related information
    *
    * @param vertices: vertices of the componnet
    * @param smoothed_v: smoothed vertices of the component
    * @param normals: normal of the component vertices
    * @param kdTree: kd-tree for neighbor query
    * @param dist_graph: [OUT] a light-weight graph with essential connection for building MST
    * @param max_length: [OUT] the distance of the longest connection each vertex involved
    * @param pre_max_length: [OUT] the maximum length of connection before connecting handles (conservative connection)
    * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
    * @param k
    * @param isEuclidean
    *
    *
    * @return None
    */
void init_graph(
    const std::vector<Point>& vertices,
    const std::vector<Point>& smoothed_v,
    const std::vector<Vec3>& normals,
    const Tree& kdTree,
    SimpGraph& dist_graph,
    std::vector<float>& max_length,
    std::vector<float>& pre_max_length,
    const float cross_conn_thresh,
    const int k,
    const bool isEuclidean)
{
    GEL_ASSERT_EQ(vertices.size(), smoothed_v.size());
    GEL_ASSERT_EQ(vertices.size(), normals.size());

    for (int i = 0; i < vertices.size(); i++) {
        dist_graph.add_node();
    }
    Util::ImmediatePool pool;
    auto neighbor_map = calculate_neighbors(pool, vertices, kdTree, k);

    NodeID i = 0;
    for (auto& vertex : vertices) {
        Vec3 this_normal = normals[i];
        auto& neighbors = neighbor_map[i];

        pre_max_length[i] = neighbors[static_cast<size_t>(k * 2. / 3.)].distance;

        // Filter out the cross-connection
        filter_out_cross_connection(neighbors, normals, i, cross_conn_thresh, isEuclidean);

        // TODO: once again, connect_nodes can be run in parallel with one protected variable
        for (const auto& neighbor : std::ranges::subrange(neighbors.begin() + 1, neighbors.end())) {
            const auto idx = neighbor.id;
            if (dist_graph.find_edge(i, idx) != AMGraph::InvalidEdgeID)
                continue;
            if (idx == i) {
                std::cerr << __func__ << ": Vertex " << idx << " connect back to its own." << std::endl;
                continue;
            }
            const Vec3 neighbor_normal = normals[idx];
            const Point neighbor_pos = vertices[idx];
            const Vec3 edge = neighbor_pos - vertex;
            const double Euclidean_dist = edge.length();
            const auto weight =
                (isEuclidean)
                    ? static_cast<float>(Euclidean_dist)
                    : static_cast<float>(cal_proj_dist(edge, this_normal, neighbor_normal));
            if (weight > max_length[i])
                max_length[i] = weight;
            if (weight > max_length[idx])
                max_length[idx] = weight;

            if (weight < 1e-8) [[unlikely]]
                std::cerr << __func__ << ": error" << std::endl;

            dist_graph.connect_nodes(i, idx, weight);
        }
        i++;
    }
}

/**
    * @brief Find the shortest path from one vertex to another in the graph
    *
    * @param mst: the graph
    * @param start: the source vertex
    * @param target: the target vertex
    * @param threshold: the step threshold, if longer than this threshold, the algorithm early stop (never used?)
    * @param path: stores the indices of vertex in the shortest path (for visualization)
    *
    * @return the number of steps of the shortest path
    */
int find_shortest_path(const RSGraph& mst, const NodeID start, const NodeID target, const int threshold,
                       std::vector<NodeID>& path)
{
    std::queue<NodeID> q;
    std::vector<int> dist(mst.no_nodes(), -1); // Distance from the start to each node
    std::vector<NodeID> pred(mst.no_nodes(), -1); // Predecessor array for path reconstruction
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
        for (const auto& v : mst.neighbors(u)) {
            if (dist[v] == -1) {
                // If the node hasn't been visited
                dist[v] = dist[u] + 1; // Increment distance
                pred[v] = u; // Record predecessor
                q.push(v);
            }
        }
    }

    return dist[target];
}

/**
* @brief weighted smoothing method using defined neighborhood with tangential distance weighted
*
* @param pool thread pool
* @param vertices: vertices of the point cloud
* @param normals: normal of the point cloud
* @param neighbors_
* @param smoothed_v: [OUT] vertices after smoothing
*
* @return None
*/
void weighted_smooth(
    Util::IExecutor& pool,
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const NeighborMap& neighbors_,
    std::vector<Point>& smoothed_v)
{
    auto lambda = [&normals, &vertices](const size_t idx, const Point& vertex, const NeighborArray& neighbors) {
        const Vec3 normal = normals[idx];

        double weight_sum = 0.;
        double amp_sum = 0.;
        double max_dist = 0.;

        std::vector<double> vertical_length;
        std::vector<double> weights;
        const std::intptr_t limit = (neighbors.size() < 192) ? static_cast<intptr_t>(neighbors.size()) : 192;
        vertical_length.reserve(limit);
        weights.reserve(limit);
        for (auto begin = neighbors.cbegin(); begin != neighbors.cbegin() + limit; ++begin) {
            const auto& neighbor = *begin;

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
                std::cerr << "error" << std::endl;
            }

            const double weight = -tangential_dist;
            if (tangential_dist > max_dist)
                max_dist = tangential_dist;

            weights.push_back(weight);
            vertical_length.push_back(vertical);
        }
        for (int i = 0; i < vertical_length.size(); i++) {
            amp_sum += vertical_length[i] * (weights[i] + max_dist);
            weight_sum += weights[i] + max_dist;
        }

        weight_sum = (weight_sum == 0.) ? 1. : weight_sum;

        amp_sum /= weight_sum;
        if (!std::isfinite(amp_sum)) [[unlikely]]
            std::cout << "error" << std::endl;
        const Vec3 move = amp_sum * normal;
        return vertex + move;
    };
    Util::Parallel::enumerate_map2(pool, vertices, neighbors_, smoothed_v, lambda);
}

auto normalize_normals(std::vector<Vec3>& normals) -> void
{
    for (auto& normal : normals) {
        normal.normalize();
    }
}

void estimate_normal_no_normals_memoized(
    Util::IExecutor& pool,
    const std::vector<Point>& vertices,
    const NeighborMap& neighbors,
    std::vector<Vec3>& normals)
{
    normals.clear();
    // Data type transfer & Cal diagonal size
    auto lambda = [&](const NeighborArray& neighbors_of_this) {
        // need id, distance and coords anyway

        std::vector<Point> neighbor_coords;
        for (const auto neighbor_id : neighbors_of_this) {
            neighbor_coords.push_back(vertices[neighbor_id.id]);
        }
        const Vec3 normal = estimateNormal(neighbor_coords);

        if (std::isnan(normal.length())) [[unlikely]] {
            std::cerr << neighbors_of_this.size() << std::endl;
            std::cerr << "error" << std::endl;
        }
        return normal;
    };
    Util::Parallel::map(pool, neighbors, normals, lambda);
}

/**
    * @brief Calculate cos angle weight for correcting normal orientation
    *
    * @param this_normal: normal of current vertex
    * @param neighbor_normal: normal of its neighbor vertex
    *
    * @return angle weight calculated
    */
double cal_angle_based_weight(const Vec3& this_normal, const Vec3& neighbor_normal)
{
    const double dot_pdt = std::abs(
        dot(this_normal, neighbor_normal) / (this_normal.length() * neighbor_normal.length()));
    const double dot_pdt_clamped = std::clamp<double>(dot_pdt, 0., 1.0);
    return 1.0 - dot_pdt_clamped;
}

/// Generic template for creating an MST
template <std::derived_from<AMGraph> TargetGraph,
    std::derived_from<AMGraph> SourceGraph,
    std::invocable<TargetGraph&, NodeID> NodeInserter, // TargetGraph -> NodeID -> ()
    std::invocable<const SourceGraph&, NodeID, NodeID> DistanceFunc, // SourceGraph -> NodeID -> NodeID -> double
    std::invocable<TargetGraph&, NodeID, NodeID> EdgeInserterFunc> // TargetGraph -> NodeID -> NodeID -> ()
TargetGraph mst(
    const SourceGraph& g,
    NodeID root,
    NodeInserter&& node_inserter,
    DistanceFunc&& distance_function,
    EdgeInserterFunc&& edge_insertion)
{
    GEL_ASSERT_NEQ(root, AMGraph::InvalidNodeID);
    TargetGraph gn;
    using QElem = std::tuple<double, NodeID, NodeID>;

    for (auto n : g.node_ids())
        node_inserter(gn, n);

    AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);

    std::priority_queue<QElem, std::vector<QElem>, std::greater<>> queue;
    for (auto neighbor : g.neighbors(root)) {
        const auto distance = distance_function(g, neighbor, root);
        queue.emplace(distance, root, neighbor);
    }

    while (!queue.empty()) {
        auto [distance, id1, id2] = queue.top();
        queue.pop();

        if (!in_tree[id2]) {
            in_tree[id2] = true;

            edge_insertion(gn, id1, id2);

            for (auto neighbor : g.neighbors(id2)) {
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
    const bool isEuclidean)
{
    GEL_ASSERT_NEQ(root, AMGraph::InvalidNodeID);
    auto gn = mst<RSGraph>(
        g, root,
        [&](RSGraph& graph, NodeID n) { graph.add_node(vertices[n], normals[n]); },
        [&](const SimpGraph& graph, NodeID n, NodeID m) {
            return CGLA::sqr_length(vertices[m] - vertices[n]);
        },
        [&](RSGraph& graph, NodeID id1, NodeID id2) {
            Vec3 edge = graph.m_vertices[id2].coords - graph.m_vertices[id1].coords;
            const auto distance =
                (isEuclidean)
                    ? edge.length()
                    : cal_proj_dist(edge, graph.m_vertices[id2].normal,
                                    graph.m_vertices[id1].normal);
            graph.add_edge(id2, id1, distance);
        }
    );
    return gn;
}

SimpGraph minimum_spanning_tree(const SimpGraph& g, NodeID root)
{
    auto gn = mst<SimpGraph>(g, root,
                             [](SimpGraph& gn, NodeID n) { gn.add_node(); },
                             [](const SimpGraph& g, NodeID n, NodeID m) { return g.get_weight(n, m); },
                             [](SimpGraph& gn, NodeID n, NodeID m) { gn.connect_nodes(n, m); });
    return gn;
}

/**
    * @brief Determine the normal orientation
    *
    * @param pool Thread pool
    * @param kdTree
    * @param in_smoothed_v
    * @param normals: [OUT] normal of the point cloud with orientation corrected
    * @param k
    *
    * @return None
    */
void correct_normal_orientation(
    Util::IExecutor& pool,
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
            const Vec3 this_normal = normals[i];

            for (const auto neighbor : std::ranges::subrange(neighbors.begin() + 1, neighbors.end())) {
                if (g_angle_temp.find_edge(i, neighbor.id) != AMGraph::InvalidEdgeID)
                    continue;
                const Vec3 neighbor_normal = normals[neighbor.id];
                const double angle_weight = cal_angle_based_weight(this_normal, neighbor_normal);

                g_angle_temp.connect_nodes(i, neighbor.id, angle_weight);
            }
        }

        return std::make_tuple(g_angle_temp, sets_temp);
    }();

    const std::vector<AMGraph::NodeSet> components_vec = connected_components(g_angle, sets);

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
            for (auto vd : mst_angle.neighbors(node_id)) {
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

/**
    * @brief Get the next neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: current outgoing branch
    *
    * @return reference to the next neighbor struct
    */
const Neighbor& successor(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.upper_bound(Neighbor(u, v, branch));
    if (iter == u.ordered_neighbors.end()) iter = u.ordered_neighbors.begin(); // Wrap around
    return (*iter); // This is honestly not good practice - ONLY modification of tree_id
}

/**
    * @brief Get last neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: current outgoing branch
    *
    * @return reference to last neighbor struct
    */
const Neighbor& predecessor(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
    if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
    return (*(std::prev(iter)));
}

void init_face_loop_label(RSGraph& g)
{
    const NodeID start_v = 0;
    NodeID last_vertex = start_v;
    auto loop_step = 0UL;
    NodeID current_vertex = g.m_vertices[start_v].ordered_neighbors.begin()->v;
    do {
        auto& next_neighbor = predecessor(g, current_vertex, last_vertex);

        next_neighbor.tree_id = g.etf.accumulate();

        last_vertex = current_vertex;
        current_vertex = next_neighbor.v;

        loop_step++;
    } while (current_vertex != g.m_vertices[start_v].ordered_neighbors.begin()->v || last_vertex != start_v);

    std::cout << "Loop step initialization finished after " + std::to_string(loop_step) + " steps." << std::endl;
}

/**
    * @brief Project a vector to a plane
    *
    * @param input: Vector to be projected
    * @param normal: normal to the plane
    *
    * @return projected Vector
    */
Vec3d projected_vector(const Vec3& input, const Vec3& normal)
{
    const Vec3d normal_normed = CGLA::normalize(normal);
    return input - CGLA::dot(input, normal_normed) * normal_normed;
}

/**
    * @brief Check if two segments are intersecting on both planes (defined by their normal) they belong
    *
    * @param mst: graph and vertex information
    * @param v1: 1st vertex of segment 1
    * @param v2: 2nd vertex of segment 1
    * @param v3: 1st vertex of segment 2
    * @param v4: 2nd vertex of segment 2
    *
    * @return if they are intersecting with each other
    */
bool is_intersecting(const RSGraph& mst, const NodeID v1, const NodeID v2, const NodeID v3, const NodeID v4)
{
    const Point p1 = mst.m_vertices[v1].coords;
    const Point p2 = mst.m_vertices[v2].coords;
    const Vec3 n1 = mst.m_vertices[v1].normal;
    const Vec3 n2 = mst.m_vertices[v2].normal;
    const Vec3 normal_12 = (n1 + n2) / 2.;

    const Point p3 = mst.m_vertices[v3].coords;
    const Point p4 = mst.m_vertices[v4].coords;
    const Vec3 n3 = mst.m_vertices[v3].normal;
    const Vec3 n4 = mst.m_vertices[v4].normal;
    const Vec3 normal_34 = (n3 + n4) / 2.;

    const auto check_intersection = [](auto p1, auto p2, auto p3, auto p4, auto normal) -> bool {
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

/**
    * @brief Geometry check for connection
    *
    * @param mst: graph and vertex information
    * @param candidate: the edge to be examed
    * @param kdTree: kd-tree for knn query
    *
    * @return if the candidate pass the check
    */
bool geometry_check(const RSGraph& mst, const TEdge& candidate, const Tree& kdTree)
{
    const NodeID v1 = candidate.first;
    const NodeID v2 = candidate.second;
    const Point p1 = mst.m_vertices[v1].coords;
    const Point p2 = mst.m_vertices[v2].coords;
    const Vec3 n1 = mst.m_vertices[v1].normal;
    const Vec3 n2 = mst.m_vertices[v2].normal;

    const Vec3 mean_normal = CGLA::normalize((n1 + n2) / 2.0);

    const Point search_center = p1 + (p2 - p1) / 2.0;
    const double radius = (p2 - p1).length() / 2.0;
    std::vector<NodeID> neighbors;
    nn_search(search_center, kdTree, radius * 3., neighbors);
    UnorderedSet<NodeID> rejection_neighbor_set;
    for (const unsigned long neighbor : neighbors) {
        if (neighbor == v1 || neighbor == v2)
            continue;
        if (CGLA::dot(mst.m_vertices[neighbor].normal, mean_normal) > std::cos(60. / 180. * M_PI))
            rejection_neighbor_set.insert(neighbor);
    }
    for (const auto neighbor : neighbors) {
        if (!rejection_neighbor_set.contains(neighbor))
            continue;

        for (const auto& rej_neighbor_neighbor : mst.m_vertices[neighbor].ordered_neighbors) {
            if (!rejection_neighbor_set.contains(rej_neighbor_neighbor.v))
                continue;

            if (is_intersecting(mst, v1, v2, neighbor, rej_neighbor_neighbor.v)) {
                return false;
            }
        }
        rejection_neighbor_set.erase(neighbor);
    }
    return true;
}

bool vanilla_check(const RSGraph& mst, const TEdge& candidate, const Tree& kdTree)
{
    const auto [this_v, neighbor] = candidate;

    // Topology check
    const auto this_v_tree = predecessor(mst, this_v, neighbor).tree_id;
    const auto neighbor_tree = predecessor(mst, neighbor, this_v).tree_id;

    if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
        return false;
    }

    return geometry_check(mst, candidate, kdTree);
}

/**
    * @brief Find the common neighbors that two vertices are sharing
    *
    * @param g: current graph
    * @param neighbor: one vertex
    * @param root: the other vertex
    * @param shared_neighbors: [OUT] common neighbors these two vertices share
    *
    * @return reference to last neighbor struct
    */
void find_common_neighbor(
    const RSGraph& g,
    const NodeID neighbor,
    const NodeID root,
    std::vector<NodeID>& shared_neighbors)
{
    shared_neighbors = g.neighbors(neighbor);
    UnorderedSet<NodeID> neighbor_neighbors(shared_neighbors.begin(), shared_neighbors.end());
    shared_neighbors = g.neighbors(root);

    // set_intersection
    std::erase_if(shared_neighbors, [&](const auto& i) {
        return !neighbor_neighbors.contains(i);
    });
}

/**
    * @brief Get the neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: the outgoing branch
    *
    * @return reference to the neighbor struct
    */
const Neighbor& get_neighbor_info(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    const auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
    return (*iter);
}

void maintain_face_loop(RSGraph& g, const NodeID source, const NodeID target)
{
    const auto this_v_tree = predecessor(g, source, target).tree_id;
    const auto neighbor_tree = predecessor(g, target, source).tree_id;

    auto [fst, snd] = g.etf.insert(this_v_tree, neighbor_tree);
    get_neighbor_info(g, source, target).tree_id = fst;
    get_neighbor_info(g, target, source).tree_id = snd;
}

bool routine_check(const RSGraph& mst, const FaceType& triangle)
{
    const Point p1 = mst.m_vertices[triangle[0]].coords;
    const Point p2 = mst.m_vertices[triangle[1]].coords;
    const Point p3 = mst.m_vertices[triangle[2]].coords;
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
        mst.m_edges[mst.find_edge(triangle[0], triangle[2])].ref_time == 2 ||
        mst.m_edges[mst.find_edge(triangle[1], triangle[2])].ref_time == 2)
        return true;
    else
        return false;
}


void add_face(RSGraph& G, const FaceType& item,
              std::vector<FaceType>& faces)
{
    const NodeID v_i = item[0];
    const NodeID v_u = item[1];
    const NodeID v_w = item[2];

    G.m_edges[G.find_edge(v_u, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_i, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_u, v_i)].ref_time += 1;
    faces.push_back(item);
}

bool register_face(RSGraph& mst, const NodeID v1, const NodeID v2, std::vector<FaceType>& faces,
                   const float edge_length)
{
    const Point p1 = mst.m_vertices[v1].coords;
    const Point p2 = mst.m_vertices[v2].coords;

    if (mst.find_edge(v1, v2) != AMGraph::InvalidEdgeID)
        return false;

    std::vector<NodeID> share_neighbors;
    find_common_neighbor(mst, v1, v2, share_neighbors);
    if (share_neighbors.empty()) {
        mst.add_edge(v1, v2, edge_length);
        maintain_face_loop(mst, v1, v2);
        return true;
    }

    const auto possible_root1 = predecessor(mst, v1, v2).v;
    const auto angle1 = cal_radians_3d(p1 - mst.m_vertices[possible_root1].coords,
                                       mst.m_vertices[possible_root1].normal,
                                       p2 - mst.m_vertices[possible_root1].coords);
    const auto possible_root2 = predecessor(mst, v2, v1).v;
    const auto angle2 = cal_radians_3d(p2 - mst.m_vertices[possible_root2].coords,
                                       mst.m_vertices[possible_root2].normal,
                                       p1 - mst.m_vertices[possible_root2].coords);

    bool isValid = true;
    std::vector<FaceType> temp;
    for (const auto v3 : share_neighbors) {
        FaceType triangle{v1, v2, v3};
        if (v3 == possible_root1 && angle1 < M_PI) {
            if (routine_check(mst, triangle)) {
                isValid = false;
                break;
            }
            if (successor(mst, v2, v1).v != v3) {
                isValid = false;
                break;
            }
            if (successor(mst, possible_root1, v2).v != v1) {
                isValid = false;
                break;
            }
            temp.emplace_back(FaceType{v1, v3, v2});
        }

        if (v3 == possible_root2 && angle2 < M_PI) {
            if (routine_check(mst, triangle)) {
                isValid = false;
                break;
            }
            if (successor(mst, v1, v2).v != v3) {
                isValid = false;
                break;
            }
            if (successor(mst, possible_root2, v1).v != v2) {
                isValid = false;
                break;
            }
            temp.emplace_back(FaceType{v1, v2, v3});
        }
    }

    if (temp.empty())
        isValid = false;

    if (isValid) {
        mst.add_edge(v1, v2, edge_length);
        maintain_face_loop(mst, v1, v2);
        for (auto& face : temp) {
            add_face(mst, face, faces);
        }
    }

    return isValid;
}

/**
    * @brief Connect handle to raise the genus number
    *
    * @param smoothed_v: smoothed vertices of the point cloud
    * @param mst: graph and vertex information
    * @param kdTree: kd-tree for knn query
    * @param connected_handle_root: [OUT] log the connected handles
    * @param k: number of kNN search
    * @param isEuclidean: if to use Euclidean distance
    * @param step_thresh: step threshold for shortest distance path early stop
    *
    * @return None
    */
void connect_handle(
    const std::vector<Point>& smoothed_v,
    Tree& kdTree,
    RSGraph& mst,
    std::vector<NodeID>& connected_handle_root,
    int k,
    int step_thresh,
    bool isEuclidean)
{
    std::vector<NodeID> imp_node;
    int num = 0;
    int edge_num = 0;
    // Collect vertices w/ an open angle larger than pi
    {
        for (int i = 0; i < mst.no_nodes(); i++) {
            const auto& neighbors = mst.m_vertices[i].ordered_neighbors;
            auto last_angle = (--neighbors.end())->angle;

            for (const auto& neighbor : neighbors) {
                auto this_angle = neighbor.angle;
                auto angle_diff = this_angle - last_angle;
                if (angle_diff < 0)
                    angle_diff += 2 * M_PI;
                if (angle_diff > M_PI)
                    imp_node.push_back(i);
                last_angle = this_angle;
            }
        }
    }

    std::vector<NodeID> connect_p;
    std::vector<NodeID> to_connect_p;
    std::vector<uint> tree_id;
    std::vector<uint> to_tree_id;

    Util::ImmediatePool pool;
    const auto neighbors_map = calculate_neighbors(pool, smoothed_v, imp_node, kdTree, k);

    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    size_t i = 0;
    for (const auto& this_v : imp_node) {
        const auto& neighbors = neighbors_map[i++];

        // Potential handle collection
        uint tree, to_tree;
        NodeID validIdx = -1;

        last_dist += (smoothed_v[this_v] - last_v).length();
        last_dist = neighbors[neighbors.size() - 1].distance;
        last_v = smoothed_v[this_v];

        for (size_t j = 0; j < neighbors.size(); j++) {
            auto& neighbor = neighbors[j];
            TEdge candidate(this_v, neighbor.id);
            if (mst.find_edge(this_v, neighbor.id) != AMGraph::InvalidEdgeID)
                continue;
            tree = mst.etf.representative((predecessor(mst, this_v, neighbor.id).tree_id));
            to_tree = mst.etf.representative(predecessor(mst, neighbor.id, this_v).tree_id);
            if (geometry_check(mst, candidate, kdTree) && tree != to_tree) {
                validIdx = j;
                break;
            }
        }
        // TODO: Check if any tree shares root, and return corresponding edges

        if (validIdx != -1) {
            connect_p.push_back(this_v);
            to_connect_p.push_back(neighbors[validIdx].id);
            tree_id.push_back(tree);
            to_tree_id.push_back(to_tree);
        }
    }

    // Select one handle
    std::unordered_map<FaceConnectionKey, std::vector<int>, FaceConnectionKeyHasher> face_connections;
    for (int i = 0; i < connect_p.size(); i++) {
        uint tree = tree_id[i];
        uint to_tree = to_tree_id[i];
        if (to_tree > tree)
            std::swap(tree, to_tree);
        auto key = FaceConnectionKey{tree, to_tree};
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
    for (auto& val : sorted_face | std::views::values) {
        const std::vector<int>& idx_vec = face_connections[val];
        if (idx_vec.size() <= 5 || (mst.exp_genus >= 0 && num >= mst.exp_genus))
            break;

        bool isFind = false;
        for (int idx : idx_vec) {
            NodeID this_v = connect_p[idx];
            Point query = mst.m_vertices[this_v].coords;
            NodeID connected_neighbor = to_connect_p[idx];
            std::vector<NodeID> path;
            const int steps = find_shortest_path(mst, this_v, connected_neighbor, step_thresh, path);
            if (steps > step_thresh) {
                isFind = true;
                const auto candidate = TEdge(this_v, connected_neighbor);
                if (geometry_check(mst, candidate, kdTree)) {
                    Vec3 edge = query - mst.m_vertices[connected_neighbor].coords;
                    const auto distance =
                        (isEuclidean)
                            ? edge.length()
                            : cal_proj_dist(edge, mst.m_vertices[this_v].normal,
                                            mst.m_vertices[connected_neighbor].normal);

                    mst.add_edge(this_v, connected_neighbor, distance);
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
    const std::vector<float>& length_thresh,
    const bool is_euclidean)
{
    const NodeID v_i = i;
    bool isFound = false;
    for (auto& neighbor : G.m_vertices[i].ordered_neighbors) {
        const NodeID v_u = neighbor.v;
        const NodeID v_w = successor(G, i, v_u).v;

        Point w_pos = G.m_vertices[v_w].coords;
        Point u_pos = G.m_vertices[v_u].coords;
        Point i_pos = G.m_vertices[v_i].coords;
        Vec3 i_normal = G.m_vertices[v_i].normal;
        Vec3 u_normal = G.m_vertices[v_u].normal;
        Vec3 w_normal = G.m_vertices[v_w].normal;
        const auto angle = cal_radians_3d(w_pos - i_pos, i_normal,
                                          u_pos - i_pos);
        const bool isLargerThanPi = angle < M_PI;
        FaceType face_vector{v_i, v_u, v_w};
        if (v_u != v_w && isLargerThanPi) {
            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                const auto score = (is_euclidean)
                                       ? (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length()
                                       : cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                                                               u_normal, w_normal);
                if (score > length_thresh[v_u] || score > length_thresh[v_w])
                    continue;
                if (score >= 0) {
                    std::pair<FaceType, float> queue_item(face_vector, score);
                    queue.push(queue_item);
                    isFound = true;
                }
            }
        }
    }

    return isFound;
}

/**
    * @brief Calculate the local surface normal (averaged direction of normals of 3 vertices in the triangle)
    *
    * @return local surface normal
    */
Vec3 triangle_mean_normal(const Vec3& normal1, const Vec3& normal2, const Vec3& normal3)
{
    Vec3 normal1_norm = normal1;
    normal1_norm.normalize();
    Vec3 normal2_norm = normal2;
    normal2_norm.normalize();
    Vec3 normal3_norm = normal3;
    normal3_norm.normalize();
    Vec3 output = normal1_norm + normal2_norm + normal3_norm;
    output.normalize();
    return output;
}

bool check_branch_validity(RSGraph& G, const NodeID root, const NodeID branch1, const NodeID branch2)
{
    constexpr auto angle_thresh = 0. / 180. * M_PI;

    const Point pos_u = G.m_vertices[branch1].coords;
    const Point pos_w = G.m_vertices[branch2].coords;
    const Vec3 normal_u = G.m_vertices[branch1].normal;
    const Vec3 normal_w = G.m_vertices[branch2].normal;

    const auto is_valid = [&root](auto this_radian, auto former, auto next) {
        bool isValid = false;
        if (next.v == root) {
            auto diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (former.v == root) {
            auto diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        return isValid;
    };

    // Check u's RS validity
    {
        const auto this_radian = cal_radians_3d(pos_w - pos_u, normal_u);
        const auto former = predecessor(G, branch1, branch2);
        const auto next = successor(G, branch1, branch2);
        if (!is_valid(this_radian, former, next))
            return false;

        // Thresh on angle
        auto diff_angle_thresh = this_radian - former.angle;
        if (diff_angle_thresh < 0)
            diff_angle_thresh += M_PI * 2.;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }


    {
        //Check w
        const auto this_radian = cal_radians_3d(pos_u - pos_w, normal_w);
        const auto former = predecessor(G, branch2, branch1);
        const auto next = successor(G, branch2, branch1);
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
    RSGraph& G,
    const std::pair<FaceType, float>& item)
{
    const NodeID v_i = item.first[0];
    const NodeID v_u = item.first[1];
    const NodeID v_w = item.first[2];
    const Point pos_i = G.m_vertices[v_i].coords;
    const Point pos_u = G.m_vertices[v_u].coords;
    const Point pos_w = G.m_vertices[v_w].coords;
    const Vec3 normal_i = G.m_vertices[v_i].normal;

    if (G.find_edge(v_u, v_w) != AMGraph::InvalidEdgeID)
        return false;

    // Non-manifold edge check
    if (G.m_edges[G.find_edge(v_i, v_u)].ref_time == 2 ||
        G.m_edges[G.find_edge(v_i, v_w)].ref_time == 2)
        return false;

    // Check this rotation system
    bool isValid = (successor(G, v_i, v_u).v == v_w);
    const auto angle = cal_radians_3d(pos_w - pos_i, normal_i, pos_u - pos_i);
    if (angle > M_PI)
        isValid = false;

    if (!isValid)
        return false;

    // Check the rotation system's validity of branch nodes
    if (!check_branch_validity(G, v_i, v_u, v_w)) {
        return false;
    }

    return true;
}

void triangulate(
    std::vector<FaceType>& faces,
    RSGraph& G,
    const bool isEuclidean,
    const std::vector<float>& length_thresh,
    const std::vector<NodeID>& connected_handle_root)
{
    UnorderedSet<NodeID> to_visit;
    FacePriorityQueue queue;

    // Init priority queue
    for (int i = 0; i < G.no_nodes(); i++) {
        to_visit.insert(i);
    }

    for (unsigned long i : connected_handle_root) {
        //bool result =
        explore(G, i, queue, length_thresh, isEuclidean);
        to_visit.erase(i);
    }

    std::cout << "Global init done :)" << std::endl;

    while (!to_visit.empty()) {
        while (!queue.empty()) {
            std::pair<FaceType, float> item = queue.top();
            queue.pop();

            if (item.second >= 0) {
                if (!check_validity(G, item))
                    continue;
            }

            // Add the edge
            const NodeID v_i = item.first[0];
            const NodeID v_u = item.first[1];
            const NodeID v_w = item.first[2];
            const Point pos_u = G.m_vertices[v_u].coords;
            const Point pos_w = G.m_vertices[v_w].coords;

            Vec3 edge = pos_u - pos_w;
            const auto distance =
                (isEuclidean)
                    ? edge.length()
                    : cal_proj_dist(edge, G.m_vertices[v_u].normal, G.m_vertices[v_w].normal);

            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                G.add_edge(v_u, v_w, distance);
                add_face(G, item.first, faces);
            } else
                continue;

            // Deal with incident triangles
            std::vector<NodeID> share_neighbors;
            find_common_neighbor(G, v_u, v_w, share_neighbors);
            for (const NodeID incident_root : share_neighbors) {
                if (incident_root == v_i)
                    continue;
                std::array face{incident_root, v_w, v_u};

                // Non-manifold edge check
                const int time1 = G.m_edges[G.find_edge(incident_root, v_u)].ref_time;
                const int time2 = G.m_edges[G.find_edge(incident_root, v_w)].ref_time;
                const int time3 = G.m_edges[G.find_edge(v_u, v_w)].ref_time;
                if (time1 == 2 || time2 == 2 || time3 == 2)
                    continue;

                add_face(G, face, faces);
            }

            to_visit.erase(v_u);
            to_visit.erase(v_w);

            // Explore and sanity check
            bool isFound = false;
            bool result = explore(G, v_u, queue, length_thresh, isEuclidean);
            isFound = isFound || result;
            result = explore(G, v_w, queue, length_thresh, isEuclidean);
            isFound = isFound || result;
        }

        if (!to_visit.empty()) {
            NodeID pick = *to_visit.begin();
            to_visit.erase(pick);
            [[maybe_unused]] bool result = explore(G, pick, queue, length_thresh, isEuclidean);
            // TODO: result is dropped?
        }
    }
}

/**
 * @brief Build minimum spanning tree (MST)
 *
 * @param out_mst: [OUT] constructed MST
 * @param g: connection information of the mst
 * @param root root node
 * @param isEuclidean: if to use Euclidean distance
 * @param vertices: coordinates of the point cloud
 * @param normals: normal of the point cloud
 *
 * @return None
 */
void build_mst(
    const SimpGraph& g,
    const NodeID root,
    RSGraph& out_mst,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    const bool isEuclidean)
{
    RSGraph temp = minimum_spanning_tree(g, root, normals, vertices, isEuclidean);

    // Fix strong ambiguous points
    if (!isEuclidean) {
        for (int i = 0; i < temp.m_edges.size(); i++) {
            const NodeID source = temp.m_edges[i].source;
            const NodeID target = temp.m_edges[i].target;
            const Vec3 normal1 = normalize(temp.m_vertices[source].normal);
            const Vec3 normal2 = normalize(temp.m_vertices[target].normal);
            const Point pos1 = temp.m_vertices[source].coords;
            const Point pos2 = temp.m_vertices[target].coords;
            if (temp.valence(source) >= 2 && temp.valence(target) >= 2)
                continue;
            const Vec3 edge = pos2 - pos1;

            Vec3 normal_sum = normal1 + normal2;
            const auto cos_angle = std::abs(dot(edge, normal_sum / normal_sum.length() / edge.length()));
            if (cos_angle > std::cos(10. / 180. * M_PI)) {
                NodeID parent;
                if (temp.valence(source) == 1) {
                    temp.m_vertices[source].normal = temp.m_vertices[target].normal;
                    parent = target;
                } else {
                    temp.m_vertices[target].normal = temp.m_vertices[source].normal;
                    parent = source;
                }

                for (const auto neighbor : temp.neighbors(parent)) {
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
    for (int i = 0; i < temp.m_vertices.size(); i++) {
        out_mst.add_node(temp.m_vertices[i].coords, temp.m_vertices[i].normal);
    }
    for (int i = 0; i < temp.m_edges.size(); i++) {
        const Vec3 edge = out_mst.m_vertices[temp.m_edges[i].source].coords -
            out_mst.m_vertices[temp.m_edges[i].target].coords;

        const double distance =
            (isEuclidean)
                ? edge.length()
                : cal_proj_dist(edge, out_mst.m_vertices[temp.m_edges[i].source].normal,
                                out_mst.m_vertices[temp.m_edges[i].target].normal);

        out_mst.add_edge(temp.m_edges[i].source, temp.m_edges[i].target, distance);
    }
}


auto estimate_normals_included_normals(
    const std::vector<Point>& vertices,
    std::vector<Vec3>& normals,
    const Distance dist) -> std::vector<Point>
{
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    Util::ImmediatePool pool;
    std::vector<Vec3> smoothed_v;
    std::vector<Point> temp;
    const int smoothing_size = std::max(static_cast<int>(static_cast<double>(vertices.size()) / 2000.), 192);

    normalize_normals(normals);
    if (dist == Distance::EUCLIDEAN) {
        smoothed_v = vertices;
    } else if (dist == Distance::NEIGHBORS) {
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
    Util::ImmediatePool pool;
    const int smoothing_size = std::max(static_cast<int>(static_cast<double>(vertices.size()) / 2000.), 192);

    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    const auto kdTree = build_kd_tree_of_indices(vertices, indices);
    auto neighbors = calculate_neighbors(pool, vertices, kdTree, smoothing_size);
    estimate_normal_no_normals_memoized(pool, vertices, neighbors, normals);

    std::vector<Point> smoothed_v;

    switch (dist) {
    case Distance::EUCLIDEAN:
        smoothed_v = vertices;
        break;
    case Distance::NEIGHBORS:
        weighted_smooth(pool, vertices, normals, neighbors, smoothed_v);
        break;
    default:
        GEL_ASSERT(false, "unreachable");
    }

    const auto temp_tree1 = build_kd_tree_of_indices(smoothed_v, indices);
    neighbors = calculate_neighbors(pool, vertices, temp_tree1, smoothing_size);
    estimate_normal_no_normals_memoized(pool, smoothed_v, neighbors, normals);

    if (dist == Distance::NEIGHBORS) {
        std::vector<Point> temp;
        temp.reserve(smoothed_v.size());
        std::swap(temp, smoothed_v);
        smoothed_v.clear();

        weighted_smooth(pool, temp, normals, neighbors, smoothed_v);

        const Tree temp_tree2 = build_kd_tree_of_indices(smoothed_v, indices);
        neighbors = calculate_neighbors(pool, vertices, temp_tree2, smoothing_size, std::move(neighbors));
        estimate_normal_no_normals_memoized(pool, smoothed_v, neighbors, normals);
    }
    GEL_ASSERT_EQ(vertices.size(), normals.size());
    return smoothed_v;
}

[[nodiscard]]
auto estimate_normals_and_smooth(
    Util::IExecutor& pool,
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

    explicit Components(
        std::vector<std::vector<Point>>&& vertices,
        std::vector<std::vector<Point>>&& smoothed_v,
        std::vector<std::vector<Vec3>>&& normals) noexcept
        : vertices{std::move(vertices)}, smoothed_v{std::move(smoothed_v)}, normals{std::move(normals)} {}

    Components() = delete;
    Components(Components& other) = delete;
};

/**
    * @brief Find the number of connected components and separate them
    *
    * @param pool: Thread pool to use
    * @param vertices: vertices of the point cloud
    * @param smoothed_v: smoothed vertices of the point cloud
    * @param normals: normal of the point cloud vertices
    * @param kdTree: kd-tree for the neighbor query
    * @param opts: theta: (cross-connection threshold) angle threshold to avoid connecting vertices on different surface
    *              r: (outlier_thresh) threshold distance(?) to remove an outlier
    *              k
    *              isEuclidean
    * @return None
    */
[[nodiscard]]
auto split_components(
    Util::IExecutor& pool,
    const Tree& kdTree,
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
    const double cross_conn_thresh = opts.theta;
    const double outlier_thresh = opts.r;
    const int k = opts.k;
    const bool isEuclidean = opts.dist == Distance::EUCLIDEAN;
    double avg_edge_length = 0;
    // TODO: can't we cache this?
    AMGraph::NodeSet sets;
    SimpGraph components;
    for (int i = 0; i < vertices.size(); i++) {
        sets.insert(components.add_node());
    }
    auto neighbor_map = calculate_neighbors(pool, vertices, kdTree, k);
    NodeID this_idx = 0;
    // Construct graph
    for (auto& neighbors : neighbor_map) {
        filter_out_cross_connection(neighbors, normals, this_idx, cross_conn_thresh, isEuclidean);

        for (auto& neighbor : neighbors) {
            NodeID idx = neighbor.id;
            double length = neighbor.distance;

            if (this_idx == idx)
                continue;

            avg_edge_length += length;

            for (int j = 0; j < k; j++) {
                if (components.find_edge(this_idx, idx) != AMGraph::InvalidEdgeID)
                    continue;
                components.connect_nodes(this_idx, idx);
            }
        }
        this_idx++;
    }
    const double thresh_r = avg_edge_length / static_cast<double>(components.no_edges()) * outlier_thresh;
    // Remove Edges Longer than the threshold
    std::vector<std::pair<NodeID, NodeID>> edge_rm_v;
    for (NodeID vertex1 = 0; vertex1 < components.no_nodes(); ++vertex1) {
        for (const auto& vertex2 : components.neighbors(vertex1)) {
            double edge_length = (vertices[vertex1] - vertices[vertex2]).length();

            if (edge_length > thresh_r) {
                edge_rm_v.emplace_back(vertex1, vertex2);
            }
        }
    }
    for (auto& [fst, snd] : edge_rm_v) {
        components.disconnect_nodes(fst, snd);
    }
    // Find Components
    std::vector<AMGraph::NodeSet> components_vec;
    components_vec = connected_components(components, sets);
    const auto num = components_vec.size();
    std::cout << "The input contains " << num << " connected components." << std::endl;
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

auto component_to_manifold(
    Util::IExecutor& pool,
    const RsROpts& opts,
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& smoothed_v) -> ::HMesh::Manifold
{
    // Insert the number_of_data_points in the tree
    const auto indices = std::ranges::iota_view(0UL, smoothed_v.size());
    Tree kdTree = build_kd_tree_of_indices(smoothed_v, indices);

    std::cout << "Init mst" << std::endl;

    // Initial Structure
    RSGraph mst;
    /*mst.init(vertices.size());*/
    std::vector<TEdge> full_edges;
    std::vector<PEdgeLength> edge_length;
    std::vector<float> connection_max_length(vertices.size(), 0.);
    std::vector<float> pre_max_length(vertices.size(), 0.);
    mst.exp_genus = opts.genus;
    {
        SimpGraph g;
        init_graph(smoothed_v, smoothed_v, normals,
                   kdTree, g, connection_max_length,
                   pre_max_length, opts.theta, opts.k, opts.dist == Distance::EUCLIDEAN);

        // Generate MST
        build_mst(g, 0, mst, normals, smoothed_v, opts.dist == Distance::EUCLIDEAN);

        // Edge arrays and sort
        for (NodeID node : g.node_ids()) {
            for (NodeID node_neighbor : g.neighbors(node)) {
                if (node < node_neighbor) {
                    Vec3 edge = smoothed_v[node] - smoothed_v[node_neighbor];
                    double len = edge.length();

                    if (opts.dist == Distance::NEIGHBORS) {
                        len = cal_proj_dist(edge, normals[node], normals[node_neighbor]);
                    }

                    if (len > pre_max_length[node] ||
                        len > pre_max_length[node_neighbor])
                        continue;
                    edge_length.emplace_back(len, full_edges.size());
                    full_edges.emplace_back(node, node_neighbor);
                }
            }
        }
        std::ranges::sort(edge_length.begin(), edge_length.end(), edge_comparator);
    }

    // Initialize face loop label
    mst.etf.reserve(6 * vertices.size());
    init_face_loop_label(mst);

    std::vector<FaceType> faces;
    // Vanilla MST imp
    // Edge connection
    for (auto& [edge_len, edge_index] : edge_length) {
        TEdge this_edge = full_edges[edge_index];

        if (mst.find_edge(this_edge.first, this_edge.second) == AMGraph::InvalidEdgeID) {
            // TODO: isAdded not checked?
            if (bool isValid = vanilla_check(mst, this_edge, kdTree)) {
                bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces,
                                             edge_len);
            }
        }
    }
    std::cout << std::endl;

    // Create handles & Triangulation
    if (opts.genus != 0) {
        std::vector<NodeID> connected_handle_root;
        connect_handle(smoothed_v, kdTree, mst, connected_handle_root, opts.k, opts.n,
                       opts.dist == Distance::EUCLIDEAN);
        triangulate(faces, mst, opts.dist == Distance::EUCLIDEAN, connection_max_length, connected_handle_root);
    }

    Manifold res;

    std::vector<size_t> flattened_face;
    for (auto& face : faces) {
        flattened_face.push_back(face[0]);
        flattened_face.push_back(face[1]);
        flattened_face.push_back(face[2]);
    }

    build_manifold(res, vertices, flattened_face, 3);

    return res;
}

auto point_cloud_to_mesh(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const RsROpts& opts) -> ::HMesh::Manifold
{
    auto opts2 = opts;
    ::HMesh::Manifold output;
    Util::RSRTimer timer;

    timer.start("Whole process");

    auto vertices_copy = vertices;
    auto normals_copy = normals;
    Util::ImmediatePool pool;
    if (!normals.empty()) {
        GEL_ASSERT_EQ(vertices.size(), normals.size(), "Vertices and normals must be the same size");
    }

    // Estimate normals & orientation & weighted smoothing
    timer.start("Estimate and smooth normals");

    std::vector<Point> in_smoothed_v = estimate_normals_and_smooth(pool, vertices_copy, normals_copy,
                                                                   opts2.dist);
    timer.end("Estimate and smooth normals");

    timer.start("Correct normal orientation");
    const auto indices = std::ranges::iota_view(0UL, in_smoothed_v.size());
    const Tree kdTree = build_kd_tree_of_indices(in_smoothed_v, indices);
    std::cout << "correct normal orientation\n";

    if (normals.empty()) {
        correct_normal_orientation(pool, kdTree, in_smoothed_v, normals_copy, opts2.k);
    }
    timer.end("Correct normal orientation");

    // Find components
    timer.start("Split components");
    std::cout << "find components\n";
    auto [component_vertices,
            component_smoothed_v,
            component_normals] =
        split_components(pool, kdTree, std::move(vertices_copy), std::move(normals_copy), std::move(in_smoothed_v),
                         opts2);
    timer.end("Split components");
    // There is no guarantee that there is more than one component, and components can
    // be highly non-uniform in terms of how many primitives they have. That means we cannot
    // rely on this loop for good parallelization opportunities.
    timer.start("Algorithm");
    for (size_t component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ...\n";

        std::vector<Point> vertices_of_this = std::move(component_vertices[component_id]);
        std::vector<Vec3> normals_of_this = std::move(component_normals[component_id]);
        std::vector<Point> smoothed_v_of_this = std::move(component_smoothed_v[component_id]);

        auto res = component_to_manifold(
            pool,
            opts2,
            vertices_of_this,
            normals_of_this,
            smoothed_v_of_this);

        output.merge(res);
    }
    timer.end("Algorithm");
    timer.end("Whole process");

    const std::string line(40, '=');
    std::cout << line << "\n\n";

    timer.show();

    return output;
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

auto point_cloud_collapse_reexpand(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const RsROpts& opts,
    const int max_iterations,
    const bool reexpand) -> Manifold
{
    Util::ImmediatePool pool;
    auto normals_copy = normals;

    auto _unused = estimate_normals_and_smooth(pool, vertices, normals_copy, opts.dist);

    const auto indices = std::ranges::iota_view(0UL, vertices.size());
    const Tree kd_tree = build_kd_tree_of_indices(vertices, indices);
    const auto neighbor_map = calculate_neighbors(pool, vertices, kd_tree, 30);

    auto collapse = collapse_points(vertices, normals_copy, max_iterations);

    const auto vertices_new = indexed_select(vertices, collapse.m_remaining);
    const auto normals_new = std::vector<Vec3>(); //indexed_select(normals_copy, collapse.m_remaining);

    auto manifold = point_cloud_to_mesh(vertices_new, normals_new, opts);

    if (reexpand)
        reexpand_points(manifold, std::move(collapse), vertices);

    return manifold;
}

auto point_cloud_normal_estimate(const std::vector<Point>& vertices,
                                 const std::vector<Vec3>& normals, const bool isEuclidean) -> NormalEstimationResult
{
    auto normals_copy = normals;
    const auto dist = isEuclidean ? Distance::EUCLIDEAN : Distance::NEIGHBORS;
    Util::ImmediatePool pool;
    const auto smoothed_v = estimate_normals_and_smooth(pool, vertices, normals_copy, dist);
    return {vertices, normals_copy, smoothed_v};
}
} // namespace GEL::HMesh::RsR
