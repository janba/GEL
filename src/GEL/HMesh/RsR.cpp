#include <sstream>
#include <map>
#include <unordered_set>
#include "RsR.h"

namespace HMesh::detail
{
static bool isGTNormal = true;
static bool isEuclidean = true;
static bool isFaceLoop = true;
static bool isDebug = false;
static bool isNoiseExperiment = false;
static int k = 30;
static double r = 20.;
static double theta = 60.;
static int step_thresh = 50;
static int exp_genus = -1;
static std::string model_path;
static std::string root_path;
static std::string model_name;
static std::string mode;
static RsR_Timer recon_timer;
static int bettiNum_1 = 0;


struct m_cmp {
    bool operator()(const std::pair<std::vector<NodeID>, float>& left,
        const std::pair<std::vector<NodeID>, float>& right) {
        return (left.second) > (right.second);
    }
};

typedef std::priority_queue<std::pair<std::vector<NodeID>, float>,
    std::vector<std::pair<std::vector<NodeID>, float>>,
    m_cmp> m_priority_queue;

typedef std::pair<float, int> m_Edge_length;
typedef std::pair<int, std::string> m_face_pair;
typedef std::pair<double, NodeID> m_neighbor_pair;

inline bool edge_comparator(const m_Edge_length& l, const m_Edge_length& r) {
    return l.first < r.first;
}
inline bool face_comparator(const m_face_pair& l, const m_face_pair& r) {
    return l.first > r.first;
}
inline bool neighbor_comparator(const m_neighbor_pair& l, const m_neighbor_pair& r) {
    return l.first > r.first;
}

void remove_duplicate_vertices(std::vector<Point>& vertices,
    std::vector<Vector>& normals, Tree& kdTree) {
    std::vector<Point> new_vertices;
    std::vector<Vector> new_normals;
    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    int this_idx = 0;
    int removed = 0;
    for (auto& vertex : vertices) {
        std::vector<NodeID> neighbors;
        std::vector<double> neighbor_distance;
        bool isInsert = true;
        if (std::isinf(vertex[0]) || std::isinf(vertex[1]) || std::isinf(vertex[2]))
            isInsert = false;
        if (std::isnan(vertex[0]) || std::isnan(vertex[1]) || std::isnan(vertex[2]))
            isInsert = false;
        if (isInsert) {
            kNN_search(vertex, kdTree, k, neighbors, neighbor_distance, last_dist, true);
            //last_v = vertex;

            for (int i = 0; i < neighbors.size(); i++) {
                NodeID idx = neighbors[i];
                double length = neighbor_distance[i];

                if (this_idx == idx)
                    continue;

                // Remove duplicate vertices
                if (length < 1e-8 && this_idx != idx) {
                    isInsert = false;
                    removed++;
                    break;
                }
            }
        }

        if (isInsert) {
            new_vertices.push_back(vertex);
            if (normals.size() == vertices.size())
                new_normals.push_back(normals[this_idx]);
        }
        this_idx++;
    }
    std::cout << removed << " duplicate vertices removed." << std::endl;
    vertices = new_vertices;
    normals = new_normals;
    return;
}
/**
    * @brief Clamp the value between the upper and lower bound
    *
    * @param value: the value to be clamped
    * @param upper_bound: upper bound
    * @param lower_bound: lower bound
    *
    * @return clamped value
    */
double clamp(double value, double upper_bound, double lower_bound) {
    double output = value;
    if (output > upper_bound)
        output = upper_bound;
    if (output < lower_bound)
        output = lower_bound;
    return output;
}

/**
    * @brief Calculate the reference vector for rotation system
    *
    * @param normal: normal direction for the target vertex
    * @param ref_vec: [OUTPUT] output the reference vector
    *
    * @return None
    */
void calculate_ref_vec(const Vector& normal, Vector& ref_vec) {
    float eps = 1e-6;
    float second = normal[1];
    if (second == 0.)
        second += eps;
    ref_vec = Vector(0, -normal[2] / second, 1);
    if (normal[2] == 1.)
        ref_vec = Vector(0., 1., 0.);
    ref_vec /= ref_vec.length();
    return;
}

/**
    * @brief Calculate the radian in the rotation system
    *
    * @param branch_vec: vector of the out-going edge
    * @param normal: normal of the root vertex
    *
    * @return radian
    */
double cal_radians_3d(const Vector& branch_vec, const Vector& normal) {
    Vector proj_vec = branch_vec - dot(normal, branch_vec) /
        normal.length() * normal;

    Vector ref_vec = Vector(0, 0, 0);
    calculate_ref_vec(normal, ref_vec);

    if (proj_vec.length() == 0.0)
        return 0.;

    Vector proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    float value = clamp(dot(proj_vec, proj_ref) / proj_vec.length() /
        proj_ref.length(), 1, -1);
    double radian = std::acos(value);
    if (dot(cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;

    if (std::isnan(radian)) {
        // Suppress per-call debug output from hot path.
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
double cal_radians_3d(const Vector& branch_vec, const Vector& normal, const Vector& ref_vec) {
    Vector proj_vec = branch_vec - dot(normal, branch_vec) /
        normal.length() * normal;
    if (std::abs(proj_vec.length()) < 1e-8)
        return 0.;

    Vector proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    float value = clamp(
        dot(proj_vec, proj_ref) / proj_vec.length() /
        proj_ref.length(), 1, -1);
    double radian = std::acos(value);
    if (dot(CGLA::cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;
    return radian;
}

//TODO: adapt kdtree
void build_KDTree(Tree& kdTree, std::vector<Point> vertices, std::vector<NodeID> indices) {
    int idx = 0;
    for (const auto& vertex : vertices) {
        kdTree.insert(vertex, indices[idx]);
        idx++;
    }
    kdTree.build();
}

/**
    * @brief k nearest neighbor search
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param num: number of nearest neighbors to be queried
    * @param neighbors: [OUT] indices of k nearest neighbors
    * @param neighbor_distance: [OUT] corresponding distance to the query point
    * @param isContain: does the query point itself count as a neighbor
    *
    * @return None
    */
// TODO: We can cache search result for every point after smoothing
void kNN_search(const Point& query, const Tree& kdTree,
    int num, std::vector<NodeID>& neighbors,
    std::vector<double>& neighbor_distance, double last_dist, bool isContain) {
    if (!isContain)
        num -= 1;
    std::vector<Record> records;
    records = kdTree.m_closest(num + 1, query, INFINITY);

    // Sort to the normal order from heap
    std::sort_heap(records.begin(), records.end());
    /*for (size_t i = records.size(); i > 0; --i) {
        std::pop_heap(records.begin(), records.begin() + i);
    }*/

    int idx = 0;
    for (auto& record : records) {
        if (idx == 0 && isContain) {
            idx++;
            continue;
        }
        neighbors.push_back(record.v);
        neighbor_distance.push_back(std::sqrt(record.d));
        idx++;
    }
    return;
}

/**
    * @brief neighbor search within a specific radius
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param radius: the radius of the search ball
    * @param neighbors: [OUT] indices of k nearest neighbors
    * @param neighbor_distance: [OUT] corresponding distance to the query point
    * @param isContain: does the query point itself count as a neighbor
    *
    * @return None
    */
void NN_search(const Point& query, const Tree& kdTree,
    double dist, std::vector<NodeID>& neighbors,
    std::vector<double>& neighbor_dist, bool isContain) {

    std::vector<Point> neighbor_coords;
    kdTree.in_sphere(query, dist, neighbor_coords, neighbors);

    std::vector<m_neighbor_pair> paired;
    int i = 0;
    for (auto& neighbor_coord : neighbor_coords) {
        neighbor_dist.push_back((neighbor_coord - query).length());
        paired.push_back({ neighbor_dist[i], neighbors[i] });
        i++;
    }
    std::sort(paired.begin(), paired.end(), neighbor_comparator);
    for (size_t i = 0; i < paired.size(); ++i) {
        if (i == 0 && isContain)
            continue;
        neighbor_dist[i] = paired[i].first;
        neighbors[i] = paired[i].second;
    }
    return;
}


/**
    * @brief Find the number of connected components and separate them
    *
    * @param vertices: vertices of the point cloud
    * @param component_vertices: [OUT] point cloud vertices in different connected components
    * @param smoothed_v: smoothed vertices of the point cloud
    * @param component_smoothed_v: [OUT] smoothed point cloud vertices in different connected components
    * @param normals: normal of the point cloud vertices
    * @param component_normals: [OUT] normal of the point cloud vertices in different components
    * @param kdTree: kd-tree for neighbor query
    * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
    * @param outlier_thresh: threshold to remove outlier
    *
    *
    * @return None
    */
float find_components(std::vector<Point>& vertices,
    std::vector<std::vector<Point>>& component_vertices,
    std::vector<Point>& smoothed_v,
    std::vector<std::vector<Point>>& component_smoothed_v,
    std::vector<Vector>& normals,
    std::vector<std::vector<Vector>>& component_normals,
    const Tree& kdTree, float cross_conn_thresh, float outlier_thresh) {

    double avg_edge_length = 0;
    AMGraph::NodeSet sets;
    SimpGraph components;
    for (int i = 0; i < vertices.size(); i++) {
        sets.insert(components.add_node());
    }

    NodeID this_idx = 0;
    std::set<NodeID> dup_remove;
    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    // Construct graph
    for (auto& vertex : smoothed_v) {
        if (dup_remove.find(this_idx) != dup_remove.end()) {
            this_idx++;
            continue;
        }

        std::vector<NodeID> neighbors;
        std::vector<double> neighbor_distance;
        last_dist += (vertex - last_v).length();
        kNN_search(vertex, kdTree, k, neighbors, neighbor_distance, last_dist, true);
        last_dist = neighbor_distance[neighbor_distance.size() - 1];
        last_v = vertex;

        // Filter out cross connection
        {
            std::vector<NodeID> temp;
            Vector this_normal = normals[this_idx];
            for (int j = 0; j < neighbors.size(); j++) {
                int idx = neighbors[j];
                Vector neighbor_normal = normals[idx];
                float cos_theta = dot(this_normal, neighbor_normal) /
                    this_normal.length() / neighbor_normal.length();
                float cos_thresh = std::cos(cross_conn_thresh / 180. * M_PI);
                if (isEuclidean)
                    cos_thresh = 0.;
                if (cos_theta >= cos_thresh) {
                    temp.push_back(idx);
                }
            }
            if (temp.size() == 0) {
                neighbors.clear();
            }
            else {
                neighbors.clear();
                neighbors = temp;
            }
        }

        for (int i = 0; i < neighbors.size(); i++) {
            NodeID idx = neighbors[i];
            double length = neighbor_distance[i];

            if (this_idx == idx)
                continue;

            // Remove duplicate vertices
            if (length < 1e-8 && this_idx != idx) {
                dup_remove.insert(idx);
            }

            avg_edge_length += length;

            for (int i = 0; i < k; i++) {
                if (components.find_edge(this_idx, idx) != AMGraph::InvalidEdgeID)
                    continue;
                components.connect_nodes(this_idx, idx);
            }
        }
        this_idx++;
    }
    std::cout << std::to_string(dup_remove.size()) << " duplicate vertices will be removed." << std::endl;
    float thresh_r = avg_edge_length / components.no_edges() * outlier_thresh;

    // Remove Edge Longer than threshold
    std::vector<NodeID> edge_rm_v_id1, edge_rm_v_id2;
    for (NodeID i = 0; i < components.no_nodes(); i++) {
        auto& edges = components.edges(i);
        for (const auto& pair : edges) {
            NodeID vertex1 = i;
            NodeID vertex2 = pair.first;
            double edge_length = (vertices[vertex1] - vertices[vertex2]).length();
            if (dup_remove.find(vertex1) != dup_remove.end()) {
                //boost::remove_edge(*eit, components);
                edge_rm_v_id1.push_back(vertex1);
                edge_rm_v_id2.push_back(vertex2);
                continue;
            }
            if (dup_remove.find(vertex2) != dup_remove.end()) {
                edge_rm_v_id1.push_back(vertex1);
                edge_rm_v_id2.push_back(vertex2);
                continue;
            }
            if (edge_length > thresh_r) {
                edge_rm_v_id1.push_back(vertex1);
                edge_rm_v_id2.push_back(vertex2);
            }
        }
    }

    for (int i = 0; i < edge_rm_v_id1.size(); i++) {
        components.disconnect_nodes(edge_rm_v_id1[i], edge_rm_v_id2[i]);
    }

    // Find Components
    std::vector<AMGraph::NodeSet> components_vec;
    
    components_vec = connected_components(components, sets);

    int num = components_vec.size();
    std::cout << "The input contains " << num << " connected components." << std::endl;

    // Valid Components and create new vectors for components
    int valid_num = 0;
    int threshold = std::min<int>(vertices.size(), 100);
    for (auto& component : components_vec) {
        if (component.size() >= threshold) {
            valid_num++;
            std::vector<Vector> this_normals;
            std::vector<Point> this_vertices;
            std::vector<Point> this_smoothed_v;
            for (const auto& element : component) {
                this_vertices.push_back(vertices[element]);
                this_smoothed_v.push_back(smoothed_v[element]);
                if (normals.size() == vertices.size())
                    this_normals.push_back(normals[element]);
            }

            component_normals.push_back(this_normals);
            component_vertices.push_back(this_vertices);
            component_smoothed_v.push_back(this_smoothed_v);
        }
    }

    std::cout << std::to_string(valid_num) << " of them will be reconstructed." << std::endl;

    components.clear();
    return thresh_r;
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
float cal_proj_dist(const Vector& edge, const Vector& this_normal, const Vector& neighbor_normal) {
    float Euclidean_dist = edge.length();
    float neighbor_normal_length = dot(edge, normalize(neighbor_normal));
    float normal_length = dot(edge, normalize(this_normal));
    float projection_dist = sqrtf((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
    projection_dist += sqrtf((Euclidean_dist * Euclidean_dist) -
        (neighbor_normal_length * neighbor_normal_length));
    projection_dist /= 2.;
    if (std::abs(dot(normalize(this_normal), normalize(neighbor_normal))) < std::cos(15. / 180. * M_PI))
        projection_dist = Euclidean_dist;
    return projection_dist;
}

/**
    * @brief initialize the graph and related information
    *
    * @param vertices: vertices of the componnet
    * @param smoothed_v: smoothed vertices of the component
    * @param normals: normal of the component vertices
    * @param kdTree: kd-tree for neighbor query
    * @param dist_graph: [OUT] a light-weight graph with essential connection for building MST
    * @param weightmap: [OUT] edge weight of dist_graph
    * @param max_length: [OUT] the distance of the longest connection each vertex involved
    * @param pre_max_length: [OUT] the maximum length of connection before connecting handles (conservative connection)
    * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
    *
    *
    * @return None
    */
void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
    const std::vector<Vector>& normals, const Tree& kdTree, SimpGraph& dist_graph,
    std::vector<float>& max_length, std::vector<float>& pre_max_length, float cross_conn_thresh) {
    for (int i = 0; i < vertices.size(); i++) {
        dist_graph.add_node();
    }

    Point last_v(0., 0., 0.);
    double last_distance = INFINITY;
    NodeID i = 0;
    for (auto& vertex : vertices) {
        Vector this_normal = normals[i];

        std::vector<NodeID> neighbors;
        std::vector<double> dists;
        last_distance += (vertex - last_v).length();
        kNN_search(smoothed_v[i], kdTree, k, neighbors, dists, last_distance, true);
        last_distance = dists[dists.size() - 1];
        last_v = vertex;
        pre_max_length[i] = dists[int(k * 2. / 3.)];

        // Filter out cross connection
        {
            std::vector<NodeID> temp;
            for (int j = 0; j < neighbors.size(); j++) {
                int idx = neighbors[j];
                Vector neighbor_normal = normals[idx];
                float cos_theta = dot(this_normal, neighbor_normal) /
                    this_normal.length() /
                    neighbor_normal.length();
                float cos_thresh = std::cos(cross_conn_thresh / 180. * M_PI);
                if (isEuclidean)
                    cos_thresh = 0.;
                if (cos_theta >= cos_thresh) {
                    temp.push_back(idx);
                }
            }
            if (temp.size() == 0) {
                // Suppress per-vertex warning spam in large reconstructions.
            }
            else {
                neighbors.clear();
                neighbors = temp;
            }
        }

        for (int j = 0; j < neighbors.size(); j++) {
            NodeID idx = neighbors[j];
            if (dist_graph.find_edge(i, idx) != AMGraph::InvalidEdgeID)
                continue;
            if (idx == i) {
                // Suppress per-edge warning spam in large reconstructions.
                continue;
            }
            Vector neighbor_normal = normals[idx];
            Point neighbor_pos = vertices[idx];
            Vector edge = neighbor_pos - vertex;
            float Euclidean_dist = edge.length();
            float weight = Euclidean_dist;
            if (!isEuclidean) {
                weight = cal_proj_dist(edge, this_normal, neighbor_normal);
            }
            if (weight > max_length[i])
                max_length[i] = weight;
            if (weight > max_length[idx])
                max_length[idx] = weight;
            if (weight < 1e-8) {
                // Suppress per-edge warning spam in large reconstructions.
            }

            dist_graph.connect_nodes(i, idx, weight);
        }
        i++;
    }
    return;
}

/**
    * @brief Find the shortest path from one vertex to another in the graph
    *
    * @param mst: the graph
    * @param start: the source vertex
    * @param target: the target vertex
    * @param threshold: the step threshold, if longer than this threshold, the algorithm early stop
    * @param path: stores the indices of vertex in the shortest path (for visualization)
    *
    * @return the number of steps of the shortest path
    */
int find_shortest_path(const RSGraph& mst, NodeID start, NodeID target, int threshold, std::vector<NodeID>& path) {
    // Init
    int shortest_path_step = -1;

    std::queue<NodeID> q;
    std::vector<int> dist(mst.no_nodes(), -1); // Distance from start to each node
    std::vector<NodeID> pred(mst.no_nodes(), -1); // Predecessor array for path reconstruction
    std::set<NodeID> visited;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        // If the target node is reached, stop early
        if (u == target) {
            break;
        }

        // If the node has already been visited, skip it
        if (visited.find(u) != visited.end()) {
            continue;
        }

        visited.insert(u);

        // Explore neighbors
        for (const auto& v : mst.neighbors(u)) {
            if (dist[v] == -1) { // If the node hasn't been visited
                dist[v] = dist[u] + 1; // Increment distance
                pred[v] = u; // Record predecessor
                q.push(v);
            }
        }
    }

    // Reconstruct path from start to target
    if (dist[target] != -1) { // Target is reachable
        for (int v = target; v != -1; v = pred[v]) {
            path.push_back(v);
        }
        reverse(path.begin(), path.end()); // Reverse to get the correct order
    }

    return dist[target];
}

/**
    * @brief weighted smoothing method using defined neighborhood with tangential distance weighted
    *
    * @param vertices: vertices of the point cloud
    * @param smoothed_v: [OUT] vertices after smoothing
    * @param normals: normal of the point cloud
    * @param kdTree: kd-tree for knn query
    * @param tr_dist: distance container
    * @param diagonal_length: the diagonal length of the point cloud
    *
    * @return None
    */
void weighted_smooth(const std::vector<Point>& vertices,
    std::vector<Point>& smoothed_v, const std::vector<Vector>& normals,
    const Tree& kdTree, float diagonal_length) {
    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    smoothed_v = std::vector<Point>(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        Point vertex = vertices[i];
        std::vector<NodeID> neighbors;
        std::vector<double> neighbor_dist;
        Vector normal = normals[i];

        int neighbor_num = 192;
        kNN_search(vertex, kdTree, neighbor_num,
            neighbors, neighbor_dist, last_dist, true);
        //last_v = vertex;

        double weight_sum = 0.;
        double amp_sum = 0.;
        double max_dist = 0.;
        int added = 0;
        std::vector<double> vertical_length;
        std::vector<double> weights;
        for (auto& neighbor : neighbors) {
            Point neighbor_pos = vertices[neighbor];
            Vector n2this = neighbor_pos - vertex;
            if (dot(normals[neighbor], normal) < std::cos(30. / 180. * M_PI)) {
                continue;
            }
            double vertical = dot(n2this, normal);
            double n_dist = (neighbor_pos - vertex).length();

            double tangential_square = n_dist * n_dist -
                vertical * vertical;
            double tangential_dist = 0.;
            if (tangential_square > 0.)
                tangential_dist = std::sqrt(tangential_square);
            if (!std::isfinite(tangential_dist)) {
                // Suppress per-neighbor debug output in tight loops.
            }
            float weight = -tangential_dist;
            if (tangential_dist > max_dist)
                max_dist = tangential_dist;

            weights.push_back(weight);
            vertical_length.push_back(vertical);
            added++;
        }
        for (int j = 0; j < vertical_length.size(); j++) {
            amp_sum += vertical_length[j] * (weights[j] + max_dist);
            weight_sum += weights[j] + max_dist;
        }

        if (weight_sum == 0.)
            weight_sum = 1.;
        amp_sum /= weight_sum;
        if (!std::isfinite(amp_sum)) {
            // Suppress per-vertex debug output in tight loops.
        }
        Vector move = amp_sum * normal;
        smoothed_v[i] = (vertex + move);
    }
}

bool examine_valid_normal(std::vector<Vector>& normals, int v_size) {
    bool isValid = normals.size() == v_size;
    if (isValid) {
        for (int i = 0; i < normals.size(); i++) {
            Vector normal = normals[i];
            if (normal.length() == 0) {
                isValid = false;
                break;
            }
            else {
                normals[i] = normalize(normals[i]);
            }
        }
    }
    return isValid;
}

void save_point_cloud_obj(const std::vector<CGLA::Vec3d>& pts,
    const std::string& filename)
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    for (const auto& p : pts)
    {
        out << "v "
            << p[0] << " "
            << p[1] << " "
            << p[2] << "\n";
    }
}

void estimate_normal(const std::vector<Point>& vertices,
    const Tree& kdTree, std::vector<Vector>& normals, std::vector<NodeID>& zero_normal_id,
    float& diagonal_length) {

    if (!isGTNormal)
        normals = std::vector<Vector>(vertices.size());
    int neighbor_num = 15;
    if (!isEuclidean)
        neighbor_num = std::min<int>(int(vertices.size() / 2000.), 192);
    // Data type transfer & Cal diagonal size
    std::vector<float> min{ INFINITY, INFINITY, INFINITY },
        max{ -INFINITY, -INFINITY, -INFINITY };

    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);

    for (int i = 0; i < vertices.size(); i++) {
        Point point = vertices[i];

        std::vector<NodeID> neighbors;
        std::vector<double> neighbor_dist;
        kNN_search(point, kdTree, neighbor_num, neighbors, neighbor_dist, last_dist, true);
        last_v = point;

        std::vector<Point> neighbor_coords;
        for (auto idx : neighbors) {
            neighbor_coords.push_back(vertices[idx]);
        }
        if (neighbor_coords.size() < neighbor_num) {
            // Suppress per-vertex debug output in tight loops.
        }
        Vector normal = estimateNormal(neighbor_coords, last_dist);
        if (std::isnan(normal.length())) {
            // Suppress per-vertex debug output in tight loops.
            save_point_cloud_obj(neighbor_coords, "error_points.obj");
        }
        normals[i] = normal;
    }
    
    for (const auto& point : vertices) {
        if (point[0] < min[0])
            min[0] = point[0];
        if (point[1] < min[1])
            min[1] = point[1];
        if (point[2] < min[2])
            min[2] = point[2];
        if (point[0] > max[0])
            max[0] = point[0];
        if (point[1] > max[1])
            max[1] = point[1];
        if (point[2] > max[2])
            max[2] = point[2];
    }
    Point min_p(min[0], min[1], min[2]), max_p(max[0], max[1], max[2]);
    diagonal_length = (max_p - min_p).length();

    return;
}

/**
    * @brief Add noise to the normal
    *
    * @param angle: the maximum angle the normal can pivot compared to the original direction
    * @param normals: normals exposed to noise
    *
    * @return None
    */
void add_normal_noise(float angle, std::vector<Vector>& normals) {
    const int seed = 3;
    std::mt19937 mt;
    mt.seed(seed);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform_distribution(0., M_PI);

    for (auto& normal : normals) {
        Vector ref1;
        calculate_ref_vec(normal, ref1);
        Vector ref2 = CGLA::cross(normal, ref1);
        ref2 /= ref2.length();
        float phi = 2 * uniform_distribution(mt);

        Vector plane_vec = std::cos(phi) * ref1 + std::sin(phi) * ref2;
        plane_vec /= plane_vec.length();
        normal = std::cos(angle) * normal + std::sin(angle) * plane_vec;
        normal /= normal.length();
    }
}

/**
    * @brief Calculate cos angle weight for correcting normal orientation
    *
    * @param this_normal: normal of current vertex
    * @param neighbor_normal: normal of its neighbor vertex
    *
    * @return angle weight calculated
    */
float cal_angle_based_weight(const Vector& this_normal, const Vector& neighbor_normal) {
    float dot_pdt = std::abs(CGLA::dot(this_normal, neighbor_normal) / (this_normal.length() * neighbor_normal.length()));
    dot_pdt = clamp(dot_pdt, 1., 0.);
    if (1. - dot_pdt < 0) {
        // Suppress per-edge debug output in tight loops.
    }
    return 1. - dot_pdt;
}

void minimum_spanning_tree(const SimpGraph& g, NodeID root,
    RSGraph& gn, std::vector<Vector>& normals, std::vector<Point>& vertices)
{
    using QElem = std::tuple<double, NodeID, NodeID>;
    if (root == AMGraph::InvalidNodeID)
        root = 0;

    for (auto n : g.node_ids())
        gn.add_node(vertices[n], normals[n]);

    Util::AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);
    in_tree[root] = true;

    std::priority_queue<QElem> Q;
    for (auto n : g.neighbors(root)) {
        auto d = CGLA::sqr_length(vertices[n] - vertices[root]);
        Q.push(std::make_tuple(-d, root, n));
    }

    while (!Q.empty()) {
        auto [d, n, m] = Q.top();
        Q.pop();

        if (!in_tree[m]) {
            in_tree[m] = true;

            Vector edge = gn.m_vertices[m].coords - gn.m_vertices[n].coords;
            float Euclidean_dist = edge.length();
            float projection_dist = cal_proj_dist(edge, gn.m_vertices[m].normal, gn.m_vertices[n].normal);
            if (std::isnan(projection_dist) || std::isnan(Euclidean_dist)) {
                // Suppress per-edge debug output in tight loops.
            }

            if (isEuclidean)
                gn.add_edge(m, n, Euclidean_dist);
            else
                gn.add_edge(m, n, projection_dist);

            //gn.connect_nodes(n, m);
            for (auto nn : g.neighbors(m)) {
                auto d_nn_m = CGLA::sqr_length(vertices[nn] - vertices[m]);
                Q.push(std::make_tuple(-d_nn_m, m, nn));
            }
        }
    }

    return;
}

void minimum_spanning_tree(const SimpGraph& g, NodeID root, SimpGraph& gn) {
    using NodeID = NodeID;
    using QElem = std::tuple<double, NodeID, NodeID>;
    if (root == AMGraph::InvalidNodeID)
        root = 0;

    for (auto n : g.node_ids())
        gn.add_node();

    Util::AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);
    in_tree[root] = true;

    std::priority_queue<QElem> Q;
    for (auto n : g.neighbors(root)) {
        auto d = g.get_weight(n, root);
        Q.push(std::make_tuple(-d, root, n));
    }

    while (!Q.empty()) {
        auto [d, n, m] = Q.top();
        Q.pop();

        if (!in_tree[m]) {
            in_tree[m] = true;
            gn.connect_nodes(n, m);
            for (auto nn : g.neighbors(m)) {
                auto d_nn_m = g.get_weight(nn, m);
                Q.push(std::make_tuple(-d_nn_m, m, nn));
            }
        }
    }

    return;
}

/**
    * @brief Determine the normal orientation
    *
    * @param G_angle graph whose edges have angle-based weight
    * @param normals: [OUT] normal of the point cloud with orientation corrected
    *
    * @return None
    */
/**
    * @brief Determine the normal orientation
    *
    * @param G_angle graph whose edges have angle-based weight
    * @param normals: [OUT] normal of the point cloud with orientation corrected
    *
    * @return None
    */
void correct_normal_orientation(std::vector<Point>& in_smoothed_v,
    Tree& kdTree, std::vector<Vector>& normals) {
    SimpGraph g_angle;
    AMGraph::NodeSet sets;
    //g_angle.init(in_pc.vertices);
    for (int i = 0; i < in_smoothed_v.size(); i++) {
        sets.insert(g_angle.add_node());
    }

    // Init angle based graph
    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    for (int i = 0; i < in_smoothed_v.size(); i++) {
        Point vertex = in_smoothed_v[i];
        Vector this_normal = normals[i];

        std::vector<NodeID> neighbors;
        std::vector<double> dists;
        last_dist += (vertex - last_v).length();
        kNN_search(vertex, kdTree, k, neighbors, dists, last_dist, false);
        last_dist = dists[dists.size() - 1];
        last_v = vertex;

        for (int j = 0; j < neighbors.size(); j++) {
            if (g_angle.find_edge(i, neighbors[j]) != AMGraph::InvalidEdgeID)
                continue;
            Vector neighbor_normal = normals[neighbors[j]];

//            float edge_weight = 0.0*cal_angle_based_weight(this_normal, neighbor_normal) + (dists[j]/last_dist);
            float edge_weight = cal_angle_based_weight(this_normal, neighbor_normal) +
                (abs(dot(vertex - in_smoothed_v[neighbors[j]], this_normal)) +
                abs(dot(vertex - in_smoothed_v[neighbors[j]], neighbor_normal)))/last_dist;
            if (i == neighbors[j] && j != 0) {
                // Suppress per-vertex debug output in tight loops.
            }
            if (edge_weight < 0) {
                // Suppress per-edge debug output in tight loops.
            }
            g_angle.connect_nodes(i, neighbors[j], edge_weight);
        }
    }

    for (int iter = 0; iter < 0; iter++) {
        auto new_normals = normals;
        for (int i = 0; i < in_smoothed_v.size(); i++) {
            Vector normal_i = normals[i];
            Vector normal = normal_i;
            for (auto j: g_angle.neighbors(i)) {
                Vector normal_j = normals[j];
                if (dot(normal_i, normal_j) < 0) {
                    normal_j = -normal_j;
                }
                double w = 1.0/(0.0000001 + g_angle.get_weight(i, j));
                normal +=  w*normal_j;
            }
            new_normals[i] = normalize(normal);
        }
        normals = new_normals;
    }   

    std::vector<AMGraph::NodeSet> components_vec;

    components_vec = connected_components(g_angle, sets);
    for (int i = 0; i < components_vec.size(); i++) {
        SimpGraph mst_angle;
        NodeID root = *components_vec[i].begin();
        minimum_spanning_tree(g_angle, root, mst_angle);

        std::vector<bool> visited_vertex(g_angle.no_nodes(), false);

        // Start from root
        std::queue<int> to_visit;
        to_visit.push(root);
        while (!to_visit.empty()) {
            NodeID node_id = to_visit.front();
            to_visit.pop();
            Vector this_normal = normals[node_id];
            auto neighbours = mst_angle.neighbors(node_id);
            for (auto vd : neighbours) {
                if (!visited_vertex[int(vd)]) {
                    to_visit.push(int(vd));
                    visited_vertex[int(vd)] = true;
                    Vector neighbor_normal = normals[int(vd)];
                    if (CGLA::dot(this_normal, neighbor_normal) < 0) {
                        normals[int(vd)] = -normals[int(vd)];
                    }
                }
            }
        }
    }
    return;
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
const Neighbor& successor(const RSGraph& g, const NodeID& root, const NodeID& branch) {
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
const Neighbor& predecessor(const RSGraph& g, const NodeID& root, const NodeID& branch) {
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.lower_bound({ u,v,static_cast<uint>(branch) });
    if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
    return (*(std::prev(iter)));
}

void init_face_loop_label(RSGraph& g) {
    NodeID start_v = 0;
    NodeID last_vertex = start_v;
    int loop_step = 0;
    NodeID current_vertex = g.m_vertices[start_v].ordered_neighbors.begin()->v;
    std::vector<int> towers;
    do {
        auto& next_neighbor = predecessor(g, current_vertex, last_vertex);

        next_neighbor.tree_id = g.etf.accumulate();

        last_vertex = current_vertex;
        current_vertex = next_neighbor.v;

        loop_step++;

    } while (current_vertex != g.m_vertices[start_v].ordered_neighbors.begin()->v || last_vertex != start_v);

    std::cout << "Loop step initialization finished after " + std::to_string(loop_step) + " steps." << std::endl;
    return;
}

/**
    * @brief Project a vector to a plane
    *
    * @param input: Vector to be projected
    * @param normal: normal to the plane
    *
    * @return projected Vector
    */
Vector projected_vector(Vector& input, Vector& normal) {
    Vector normal_normed = normal;
    normal_normed.normalize();
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
bool isIntersecting(RSGraph& mst, NodeID v1, NodeID v2, NodeID v3, NodeID v4) {
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;
    Vector n1 = mst.m_vertices[v1].normal;
    Vector n2 = mst.m_vertices[v2].normal;
    Point midpoint_12 = p1 + (p2 - p1) / 2.;
    Vector normal_12 = (n1 + n2) / 2.;

    Point p3 = mst.m_vertices[v3].coords;
    Point p4 = mst.m_vertices[v4].coords;
    Vector n3 = mst.m_vertices[v3].normal;
    Vector n4 = mst.m_vertices[v4].normal;
    Point midpoint_34 = p3 + (p4 - p3) / 2.;
    Vector normal_34 = (n3 + n4) / 2.;


    // On the plane of edge 12
    {
        bool isIntersecting = true;
        Vector edge1 = p1 - midpoint_12;
        Vector edge2 = p3 - midpoint_12;
        Vector edge3 = p4 - midpoint_12;
        Vector proj_edge1 = projected_vector(edge1, normal_12);
        Vector proj_edge2 = projected_vector(edge2, normal_12);
        Vector proj_edge3 = projected_vector(edge3, normal_12);
        Vector pro1 = CGLA::cross(proj_edge2, proj_edge1);
        Vector pro2 = CGLA::cross(proj_edge3, proj_edge1);
        if (CGLA::dot(pro1, pro2) > 0)
            isIntersecting = false;
        if (isIntersecting) {
            edge1 = p3 - midpoint_34;
            edge2 = p1 - midpoint_34;
            edge3 = p2 - midpoint_34;
            proj_edge1 = projected_vector(edge1, normal_12);
            proj_edge2 = projected_vector(edge2, normal_12);
            proj_edge3 = projected_vector(edge3, normal_12);
            pro1 = CGLA::cross(proj_edge2, proj_edge1);
            pro2 = CGLA::cross(proj_edge3, proj_edge1);
            if (CGLA::dot(pro1, pro2) > 0)
                isIntersecting = false;
        }
        if (isIntersecting)
            return true;
    }

    // On the plane of edge 34
    {
        bool isIntersecting = true;
        Vector edge1 = p1 - midpoint_12;
        Vector edge2 = p3 - midpoint_12;
        Vector edge3 = p4 - midpoint_12;
        Vector proj_edge1 = projected_vector(edge1, normal_34);
        Vector proj_edge2 = projected_vector(edge2, normal_34);
        Vector proj_edge3 = projected_vector(edge3, normal_34);
        Vector pro1 = CGLA::cross(proj_edge2, proj_edge1);
        Vector pro2 = CGLA::cross(proj_edge3, proj_edge1);
        if (CGLA::dot(pro1, pro2) > 0)
            isIntersecting = false;
        if (isIntersecting) {
            edge1 = p3 - midpoint_34;
            edge2 = p1 - midpoint_34;
            edge3 = p2 - midpoint_34;
            proj_edge1 = projected_vector(edge1, normal_34);
            proj_edge2 = projected_vector(edge2, normal_34);
            proj_edge3 = projected_vector(edge3, normal_34);
            pro1 = CGLA::cross(proj_edge2, proj_edge1);
            pro2 = CGLA::cross(proj_edge3, proj_edge1);
            if (CGLA::dot(pro1, pro2) > 0)
                isIntersecting = false;
        }
        if (isIntersecting)
            return true;
    }
    return false;
}

/**
    * @brief Geometry check for connection
    *
    * @param mst: graph and vertex information
    * @param candidate: the edge to be examed
    * @param kdTree: kd-tree for knn query
    * @param tr_dist: distance container
    *
    * @return if the candidate pass the check
    */
bool geometry_check(RSGraph& mst, m_Edge& candidate, Tree& kdTree) {
    NodeID v1 = candidate.first;
    NodeID v2 = candidate.second;
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;
    Vector n1 = mst.m_vertices[v1].normal;
    Vector n2 = mst.m_vertices[v2].normal;

    Vector mean_normal = (n1 + n2) / 2.;
    mean_normal.normalize();

    Point search_center = p1 + (p2 - p1) / 2.;
    float radius = (p2 - p1).length() / 2.;
    std::vector<NodeID> neighbors;
    std::vector<double> distance;
    NN_search(search_center, kdTree, float(radius * 3.), neighbors, distance, false);

    float query_radian1 = cal_radians_3d(p1 - search_center, mean_normal);
    float query_radian2 = cal_radians_3d(p2 - search_center, mean_normal);
    std::set<int> rejection_neighbor_set;
    for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] == v1 || neighbors[i] == v2)
            continue;
        if (CGLA::dot(mst.m_vertices[neighbors[i]].normal, mean_normal) > std::cos(60. / 180. * M_PI))
            rejection_neighbor_set.insert(neighbors[i]);
    }
    for (int i = 0; i < neighbors.size(); i++) {
        if (rejection_neighbor_set.find(neighbors[i]) ==
            rejection_neighbor_set.end())
            continue;
        //	continue;
        NodeID rejection_neighbor = neighbors[i];
        Point rej_neighbor_pos = mst.m_vertices[rejection_neighbor].coords;
        float min_radian, max_radian;

        for (auto& rej_neighbor_neighbor : mst.m_vertices[rejection_neighbor].ordered_neighbors) {
            if (rejection_neighbor_set.find(rej_neighbor_neighbor.v) ==
                rejection_neighbor_set.end())
                continue;

            bool result = isIntersecting(mst, v1, v2, rejection_neighbor, rej_neighbor_neighbor.v);
            if (result) {
                return false;
            }
        }
        rejection_neighbor_set.erase(neighbors[i]);
    }
    return true;
}

bool Vanilla_check(RSGraph& mst, m_Edge& candidate, Tree& kdTree) {
    NodeID neighbor = candidate.second;
    NodeID this_v = candidate.first;
    Vector this_normal = mst.m_vertices[this_v].normal;
    this_normal.normalize();
    Vector neighbor_normal = mst.m_vertices[neighbor].normal;
    neighbor_normal.normalize();

    auto this_v_tree = predecessor(mst, this_v, neighbor).tree_id;
    auto neighbor_tree = predecessor(mst, neighbor, this_v).tree_id;
    if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
        return false;
    }

    return geometry_check(mst, candidate, kdTree);
}

/**
    * @brief Find the common neighbors that two vertices are sharing
    *
    * @param neighbor: one vertex
    * @param root: the other vertex
    * @param share_neighbor: [OUT] common neighbors these two vertices share
    * @param mst: tcurrent graph
    *
    * @return reference to last neighbor struct
    */
void find_common_neighbor(NodeID neighbor, NodeID root,
    std::vector<NodeID>& share_neighbor, RSGraph& g) {
    auto adj = g.neighbors(neighbor);
    std::set<NodeID> neighbor_neighbor(adj.begin(), adj.end());
    adj = g.neighbors(root);
    std::set<NodeID> vertex_neighbor(adj.begin(), adj.end());
    std::set_intersection(vertex_neighbor.begin(), vertex_neighbor.end(),
        neighbor_neighbor.begin(), neighbor_neighbor.end(),
        std::back_inserter(share_neighbor));
    return;
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
const Neighbor& get_neighbor_info(const RSGraph& g, const NodeID& root, const NodeID& branch) {
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.lower_bound({ u,v,static_cast<uint>(branch) });
    return (*iter);
}

void maintain_face_loop(RSGraph& g, const NodeID source, const NodeID target) {
    auto this_v_tree = predecessor(g, source, target).tree_id;
    auto neighbor_tree = predecessor(g, target, source).tree_id;

    auto result = g.etf.insert(this_v_tree, neighbor_tree);
    auto u = result.first;
    auto v = result.second;
    get_neighbor_info(g, source, target).tree_id = u;
    get_neighbor_info(g, target, source).tree_id = v;

    bettiNum_1++;
    return;
}

bool routine_check(RSGraph& mst, std::vector<NodeID>& triangle) {

    NodeID v1 = triangle[0];
    NodeID v2 = triangle[1];
    NodeID v3 = triangle[2];
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;
    Point p3 = mst.m_vertices[v3].coords;
    Vector n1 = mst.m_vertices[v1].normal;
    Vector n2 = mst.m_vertices[v2].normal;
    Vector n3 = mst.m_vertices[v3].normal;

    //bool isValid = (n1 * face_normal * (n2 * face_normal) < 0 ||
    //		n1 * face_normal * (n3 * face_normal) < 0);

    {
        float len_ui = (p1 - p2).length();
        float len_wi = (p3 - p2).length();
        float len_uw = (p1 - p3).length();

        float max_value = std::acos(clamp(
            dot((p3 - p2), (p1 - p2)) / len_ui /
            len_wi, 1, -1));
        float radian = std::acos(clamp(
            dot((p2 - p1), (p3 - p1)) / len_ui /
            len_uw, 1, -1));
        if (radian > max_value)
            max_value = radian;
        radian = std::acos(clamp(
            dot((p1 - p3), (p2 - p3)) / len_uw /
            len_wi, 1, -1));
        if (radian > max_value)
            max_value = radian;
        if (max_value > 175. / 180. * M_PI)
            return true;
    }


    if (mst.m_edges[mst.find_edge(v1, v3)].ref_time == 2 ||
        mst.m_edges[mst.find_edge(v2, v3)].ref_time == 2)
        return true;

    return false;
}


void add_face(RSGraph& G, std::vector<NodeID>& item,
    std::vector<std::vector<NodeID>>& faces) {
    NodeID v_i = item[0];
    NodeID v_u = item[1];
    NodeID v_w = item[2];

    // Maintain face exist for detecting holes
    {
        get_neighbor_info(G, v_i, v_u).faceExist = true;
        get_neighbor_info(G, v_u, v_w).faceExist = true;
        get_neighbor_info(G, v_w, v_i).faceExist = true;
    }

    G.m_edges[G.find_edge(v_u, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_i, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_u, v_i)].ref_time += 1;
    faces.push_back(item);
    bettiNum_1--;
    return;
}

bool register_face(RSGraph& mst, NodeID v1, NodeID v2, std::vector<std::vector<NodeID>>& faces,
    Tree& KDTree, float edge_length) {

    Vector v1_n = mst.m_vertices[v1].normal;
    Vector v2_n = mst.m_vertices[v2].normal;
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;

    if (mst.find_edge(v1, v2) != AMGraph::InvalidEdgeID)
        return false;

    std::vector<NodeID> share_neighbors;
    find_common_neighbor(v1, v2, share_neighbors, mst);
    if (share_neighbors.size() == 0) {
        mst.add_edge(v1, v2, edge_length);
        maintain_face_loop(mst, v1, v2);
        return true;
    }

    auto possible_root1 = predecessor(mst, v1, v2).v;
    float angle1 = cal_radians_3d(p1 - mst.m_vertices[possible_root1].coords,
        mst.m_vertices[possible_root1].normal, p2 - mst.m_vertices[possible_root1].coords);
    auto possible_root2 = predecessor(mst, v2, v1).v;
    float angle2 = cal_radians_3d(p2 - mst.m_vertices[possible_root2].coords,
        mst.m_vertices[possible_root2].normal, p1 - mst.m_vertices[possible_root2].coords);

    bool isValid = true;
    std::vector<std::vector<NodeID>> temp;
    for (auto v3 : share_neighbors) {
        std::vector<NodeID> triangle{ v1,v2,v3 };
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
            /*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
                isValid = false;
                break;
            }*/
            temp.push_back(std::vector<NodeID>{v1, v3, v2});
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
            /*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
                isValid = false;
                break;
            }*/
            temp.push_back(std::vector<NodeID>{v1, v2, v3});
        }
    }

    if (temp.size() == 0)
        isValid = false;

    if (isValid) {
        AMGraph::EdgeID added_edge = mst.add_edge(v1, v2,
            edge_length);
        maintain_face_loop(mst, v1, v2);
        for (auto& face : temp) {
            ;
            add_face(mst, face, faces);
        }
    }

    return isValid;
}

void export_obj(std::vector<Point>& in_vertices, RSGraph& g, std::string out_path,
    std::vector<std::vector<NodeID>>& faces) {
    std::ofstream file(out_path);
    // Write vertices
    file << "# List of geometric vertices" << std::endl;
    for (int i = 0; i < in_vertices.size(); i++) {
        Point this_coords = in_vertices[i];
        file << "v " << std::to_string(this_coords[0])
            << " " << std::to_string(this_coords[1])
            << " " << std::to_string(this_coords[2]) << std::endl;
    }
    // Write vertex normal
    file << std::endl;
    file << "# List of vertex normals" << std::endl;
    for (int i = 0; i < g.no_nodes(); i++) {
        Vector this_normal = g.m_vertices[i].normal;
        file << "vn " << std::to_string(this_normal[0])
            << " " << std::to_string(this_normal[1])
            << " " << std::to_string(this_normal[2]) << std::endl;
    }

    // Write faces
    file << std::endl;
    file << "# Polygonal face element" << std::endl;
    for (auto face : faces) {
        file << "f " << std::to_string(face[0] + 1)
            << " " << std::to_string(face[1] + 1)
            << " " << std::to_string(face[2] + 1) << std::endl;
    }
    file.close();
    return;
}

/**
    * @brief Connect handle to raise the genus number
    *
    * @param smoothed_v: smoothed vertices of the point cloud
    * @param mst: graph and vertex information
    * @param kdTree: kd-tree for knn query
    * @param tr_dist: distance container
    * @param connected_handle_root: [OUT] log the connected handles
    * @param betti: [OUT] log betti number changes
    * @param k: number of kNN search
    * @param isEuclidean: if to use Euclidean distance
    * @param step_thresh: step threshold for shortest distance path early stop
    *
    * @return None
    */
void connect_handle(const std::vector<Point>& smoothed_v, Tree& KDTree,
    RSGraph& mst, std::vector<NodeID>& connected_handle_root,
    std::vector<int>& betti) {

    std::vector<NodeID> imp_node;
    int num = 0;
    int edge_num = 0;
    // Collect vertices w/ an open angle larger than pi
    {
        for (int i = 0; i < mst.no_nodes(); i++) {
            std::set<Neighbor>& neighbors = mst.m_vertices[i].ordered_neighbors;
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
        }
    }

    std::vector<NodeID> connect_p;
    std::vector<NodeID> to_connect_p;
    std::vector<uint> tree_id;
    std::vector<uint> to_tree_id;

    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    for (auto& this_v : imp_node) {
        Point query = mst.m_vertices[this_v].coords;
        Vector query_normal = mst.m_vertices[this_v].normal;
        std::vector<NodeID> neighbors;
        std::vector<double> dists;

        // Potential handle collection
        uint tree, to_tree;
        int validIdx = -1;

        last_dist += (smoothed_v[int(this_v)] - last_v).length();
        kNN_search(smoothed_v[int(this_v)], KDTree, k, neighbors, dists, last_dist, true);
        last_dist = dists[dists.size() - 1];
        last_v = smoothed_v[int(this_v)];

        for (int i = 0; i < neighbors.size(); i++) {
            int neighbor = neighbors[i];
            m_Edge candidate(this_v, neighbor);
            if (mst.find_edge(this_v, neighbor) != AMGraph::InvalidEdgeID)
                continue;
            tree = mst.etf.representative((predecessor(mst, this_v, neighbor).tree_id));
            to_tree = mst.etf.representative(predecessor(mst, neighbor, this_v).tree_id);
            if (geometry_check(mst, candidate, KDTree) && tree != to_tree) {
                validIdx = i;
                break;
            }
        }
        // TODO: Check if any tree shares root, and return corresponding edges

        if (validIdx != -1) {
            connect_p.push_back(this_v);
            to_connect_p.push_back(neighbors[validIdx]);
            tree_id.push_back(tree);
            to_tree_id.push_back(to_tree);
        }
    }

    // Select one handle
    std::map<std::string, std::vector<int>> face_connections;
    for (int i = 0; i < connect_p.size(); i++) {
        uint tree = tree_id[i];
        uint to_tree = to_tree_id[i];
        if (to_tree > tree)
            std::swap(tree, to_tree);
        std::string key = std::to_string(tree) + "+" + std::to_string(to_tree);
        if (face_connections.find(key) == face_connections.end())
            face_connections[key] = std::vector<int>{ i };
        else {
            face_connections[key].push_back(i);
        }
    }

    // Sort
    std::vector<m_face_pair> sorted_face;
    for (auto key = face_connections.begin(); key != face_connections.end(); key++) {
        int length = face_connections[key->first].size();
        sorted_face.push_back(m_face_pair(length, key->first));
    }
    std::sort(sorted_face.begin(), sorted_face.end(), face_comparator);
    for (int i = 0; i < sorted_face.size(); i++) {
        std::string key = sorted_face[i].second;
        std::vector<int> idx_vec = face_connections[key];
        if (idx_vec.size() <= 5)
            break;
        if (mst.exp_genus >= 0 && num >= mst.exp_genus)
            break;
        Point query;
        NodeID connected_neighbor, this_v;
        Edge added_edge;
        bool isFind = false;
        for (int idx : idx_vec) {
            this_v = connect_p[idx];
            query = mst.m_vertices[this_v].coords;
            connected_neighbor = to_connect_p[idx];
            std::vector<NodeID> path;
            int steps = find_shortest_path(mst, this_v, connected_neighbor, step_thresh, path);
            if (steps > step_thresh) {
                isFind = true;
                m_Edge candidate(this_v, connected_neighbor);
                if (geometry_check(mst, candidate, KDTree)) {
                    Vector edge = query - mst.m_vertices[connected_neighbor].coords;
                    float Euclidean_dist = edge.length();
                    float projection_dist = cal_proj_dist(edge, mst.m_vertices[this_v].normal,
                        mst.m_vertices[connected_neighbor].normal);

                    if (isEuclidean) {
                        mst.add_edge(this_v, connected_neighbor, Euclidean_dist);
                        connected_handle_root.push_back(this_v);
                        connected_handle_root.push_back(connected_neighbor);
                    }
                    else {
                        mst.add_edge(this_v, connected_neighbor, projection_dist);
                        connected_handle_root.push_back(this_v);
                        connected_handle_root.push_back(connected_neighbor);
                    }

                    bettiNum_1++;
                    edge_num++;
                    betti.push_back(bettiNum_1);
                    //fs::path out_edge_path("C:/Projects_output/letters/edge_" + std::to_string(edge_num) + ".obj");
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
    return;
}

void export_edges(RSGraph& g, std::vector<NodeID>& roots, std::string out_path) {
    std::ofstream file(out_path);
    // Write vertices
    file << "# List of geometric vertices" << std::endl;
    for (int i = 0; i < roots.size(); i++) {
        Point this_coords = g.m_vertices[i].coords;
        file << "v " << std::to_string(this_coords[0])
            << " " << std::to_string(this_coords[1])
            << " " << std::to_string(this_coords[2]) << std::endl;
    }

    // Write lines
    file << std::endl;
    file << "# Line element" << std::endl;
    for (int i = 0; i < roots.size(); i += 2)
        file << "l " << i + 1
        << " " << i + 2 << std::endl;
    file.close();
    return;
}

bool explore(RSGraph& G, int i, m_priority_queue& queue,
    std::unordered_set<std::string>& faces_in_queue, float avg_edge_length
    , std::vector<float>& length_thresh) {
    NodeID v_i = i;
    bool isFound = false;
    for (auto& neighbor : G.m_vertices[i].ordered_neighbors)
    {
        NodeID v_u = neighbor.v;
        NodeID v_w = successor(G, i, v_u).v;

        Point w_pos = G.m_vertices[v_w].coords;
        Point u_pos = G.m_vertices[v_u].coords;
        Point i_pos = G.m_vertices[v_i].coords;
        Vector i_normal = G.m_vertices[v_i].normal;
        Vector u_normal = G.m_vertices[v_u].normal;
        Vector w_normal = G.m_vertices[v_w].normal;
        float angle = cal_radians_3d(w_pos - i_pos, i_normal,
            u_pos - i_pos);
        bool isLargerThanPi = angle < M_PI;
        std::vector<NodeID> face_vector{ v_i, v_u, v_w };
        if (v_u != v_w && isLargerThanPi) {
            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                float score = (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length();
                if (!G.isEuclidean) {
                    score = cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                        u_normal, w_normal);
                }
                if (score > length_thresh[v_u] || score > length_thresh[v_w])
                    continue;
                if (score >= 0) {
                    std::pair<std::vector<NodeID>, float> queue_item(face_vector, score);
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
Vector triangle_mean_normal(const Vector& normal1, const Vector& normal2, const Vector& normal3) {
    Vector normal1_norm = normal1;
    normal1_norm.normalize();
    Vector normal2_norm = normal2;
    normal2_norm.normalize();
    Vector normal3_norm = normal3;
    normal3_norm.normalize();
    Vector output = normal1_norm + normal2_norm + normal3_norm;
    output.normalize();
    return output;
}

bool check_branch_validity(RSGraph& G, NodeID root, NodeID branch1, NodeID branch2) {
    Point pos_i = G.m_vertices[root].coords;
    Point pos_u = G.m_vertices[branch1].coords;
    Point pos_w = G.m_vertices[branch2].coords;
    Vector normal_i = G.m_vertices[root].normal;
    Vector normal_u = G.m_vertices[branch1].normal;
    Vector normal_w = G.m_vertices[branch2].normal;

    // Option 1
    std::vector<Point> triangle_pos{ pos_i, pos_u, pos_w };
    Vector face_normal = triangle_mean_normal(normal_i, normal_u, normal_w);

    float angle_thresh = 0. / 180. * M_PI;

    // Check u's RS validity
    float this_radian = cal_radians_3d(pos_w - pos_u, normal_u);
    auto former = predecessor(G, branch1, branch2);
    auto next = successor(G, branch1, branch2);
    if (G.isFinal) {
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    }
    else {
        float diff = next.angle - this_radian;
        if (diff < 0)
            diff += 2 * M_PI;
        if (next.v != root || diff > M_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = this_radian - former.angle;
        if (diff_angle_thresh < 0)
            diff_angle_thresh += M_PI * 2.;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }

    //Check w
    this_radian = cal_radians_3d(pos_u - pos_w, normal_w);
    former = predecessor(G, branch2, branch1);
    next = successor(G, branch2, branch1);
    if (G.isFinal) {
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    }
    else {
        float diff = -former.angle + this_radian;
        if (diff < 0)
            diff += 2 * M_PI;
        if (former.v != root || diff > M_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = -this_radian + next.angle + 2. * M_PI;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }
    return true;
}


bool check_validity(RSGraph& G, std::pair<std::vector<NodeID>, float>& item,
    const Tree& KDTree, bool isFaceloop, bool isFinalize) {
    int i = item.first[0];
    int u = item.first[1];
    int w = item.first[2];
    NodeID v_i = i;
    NodeID v_u = u;
    NodeID v_w = w;
    Point pos_i = G.m_vertices[v_i].coords;
    Point pos_u = G.m_vertices[v_u].coords;
    Point pos_w = G.m_vertices[v_w].coords;
    Vector normal_i = G.m_vertices[v_i].normal;
    Vector normal_u = G.m_vertices[v_u].normal;
    Vector normal_w = G.m_vertices[v_w].normal;


    if (G.find_edge(v_u, v_w) != AMGraph::InvalidEdgeID)
        return false;

    // Non-manifold edge check
    if (G.m_edges[G.find_edge(v_i, v_u)].ref_time == 2 ||
        G.m_edges[G.find_edge(v_i, v_w)].ref_time == 2)
        return false;

    // Check this rotation system
    bool isValid = (successor(G, v_i, v_u).v == v_w);
    float angle = cal_radians_3d(pos_w - pos_i, normal_i, pos_u - pos_i);
    if (angle > M_PI)
        isValid = false;

    if (!isValid)
        return false;

    if (!isFinalize) {
        // Check rotation system's validity of branch nodes
        if (!check_branch_validity(G, v_i, v_u, v_w)) {
            return false;
        }
    }

    // Check face overlap
    /*if (check_face_overlap(G, item.first, KDTree, tr_dist))
        return false;*/

    return true;
}

void checkAndForce(NodeID v_u, NodeID v_w, RSGraph& G, m_priority_queue& queue,
    std::vector<float>& length_thresh) {
    std::vector<NodeID> check_v{ v_u,v_w };
    for (int i = 0; i < check_v.size(); i++) {
        NodeID v_i = check_v[i];
        std::vector<int> connection_count(G.m_vertices[v_i].ordered_neighbors.size(), 0);
        std::vector<NodeID> neighbor_ids;
        std::set<NodeID> neighbor_ids_set;

        for (auto& neighbor : G.m_vertices[v_i].ordered_neighbors) {
            neighbor_ids_set.insert(neighbor.v);
            neighbor_ids.push_back(neighbor.v);
        }

        int idx = 0;
        for (auto& neighbor : G.m_vertices[v_i].ordered_neighbors) {
            for (auto& neighbor_neighbor : G.m_vertices[neighbor.v].ordered_neighbors) {
                if (neighbor_ids_set.find(neighbor_neighbor.v) != neighbor_ids_set.end()) {
                    connection_count[idx] += 1;
                }
            }
            idx++;
        }

        std::vector<NodeID> not_full;
        for (int i = 0; i < connection_count.size(); i++) {
            if (connection_count[i] < 2) {
                not_full.push_back(i);
            }
        }

        if (not_full.size() == 2) {
            NodeID v_u = neighbor_ids[not_full[0]];
            NodeID v_w = neighbor_ids[not_full[1]];
            Vector u_normal = G.m_vertices[v_u].normal;
            Vector w_normal = G.m_vertices[v_w].normal;
            float score = (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length();
            if (!G.isEuclidean) {
                score = cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                    u_normal, w_normal);
            }
            if (score > length_thresh[v_u] || score > length_thresh[v_w])
                return;
            std::vector<NodeID> face_vector{ v_i, neighbor_ids[not_full[0]], neighbor_ids[not_full[1]] };
            std::pair<std::vector<NodeID>, float> queue_item(face_vector, -1);
            queue.push(queue_item);
            break;
        }
    }
    return;
}

void triangulate(std::vector<std::vector<NodeID>>& faces, RSGraph& G, SimpGraph& g_org,
    const Tree& KDTree, bool isFaceLoop, bool isEuclidean
    , std::vector<float>& length_thresh, std::vector<NodeID>& connected_handle_root,
    std::vector<int>& betti, bool isFinalize = false) {

    std::unordered_set<std::string> faces_in_queue;
    std::unordered_set<NodeID> to_visit;

    m_priority_queue queue;

    float avg_edge_length = G.total_edge_length / G.no_edges();

    // Init priority queue
    for (int i = 0; i < G.no_nodes(); i++) {
        to_visit.insert(i);
    }

    for (int i = 0; i < connected_handle_root.size(); i++) {
        bool result = explore(G, connected_handle_root[i], queue, faces_in_queue, avg_edge_length, length_thresh);
        to_visit.erase(connected_handle_root[i]);
    }

    std::cout << "Global init done :)" << std::endl;

    int loop_time = 0;

    while (to_visit.size() > 0) {
        while (!queue.empty()) {
            loop_time += 1;
            std::pair<std::vector<NodeID>, float> item = queue.top();
            queue.pop();

            if (item.second >= 0) {
                // Validity check
                bool isValid = check_validity(G, item, KDTree, isFaceLoop, isFinalize);
                if (!isValid)
                    continue;
            }

            // Add the edge
            NodeID v_i = item.first[0];
            NodeID v_u = item.first[1];
            NodeID v_w = item.first[2];
            Point pos_i = G.m_vertices[v_i].coords;
            Point pos_u = G.m_vertices[v_u].coords;
            Point pos_w = G.m_vertices[v_w].coords;
            Vector normal_i = G.m_vertices[v_i].normal;
            Vector normal_u = G.m_vertices[v_u].normal;
            Vector normal_w = G.m_vertices[v_w].normal;

            float dist = (pos_u - pos_w).length();
            Vector edge = pos_u - pos_w;
            float Euclidean_dist = (edge).length();
            float projection_dist = cal_proj_dist(edge,
                G.m_vertices[v_u].normal, G.m_vertices[v_w].normal);

            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                Edge added_edge;
                if (isEuclidean)
                    G.add_edge(v_u, v_w, Euclidean_dist);
                else
                    G.add_edge(v_u, v_w, projection_dist);

                avg_edge_length = G.total_edge_length / G.no_edges();

                bettiNum_1++;

                add_face(G, item.first, faces);

                betti.push_back(bettiNum_1);
            }
            else
                continue;

            // Deal with incident triangles
            {
                std::vector<NodeID> share_neighbors;
                find_common_neighbor(v_u, v_w, share_neighbors, G);
                for (int idx = 0; idx < share_neighbors.size(); idx++) {
                    NodeID incident_root = share_neighbors[idx];
                    if (incident_root == v_i)
                        continue;
                    std::vector<NodeID> face{ incident_root,v_w,v_u };

                    // Non-manifold edge check
                    int time1 = G.m_edges[G.find_edge(incident_root, v_u)].ref_time;
                    int time2 = G.m_edges[G.find_edge(incident_root, v_w)].ref_time;
                    int time3 = G.m_edges[G.find_edge(v_u, v_w)].ref_time;
                    if (time1 == 2 || time2 == 2 || time3 == 2)
                        continue;

                    add_face(G, face, faces);
                }
            }

            to_visit.erase(v_u);
            to_visit.erase(v_w);

            // Explore and sanity check
            bool isFound = false;
            bool result = explore(G, v_u, queue, faces_in_queue, avg_edge_length, length_thresh);
            isFound = isFound || result;
            result = explore(G, v_w, queue, faces_in_queue, avg_edge_length, length_thresh);
            isFound = isFound || result;

            if (isFinalize) {
                if ((!isFound)) {
                    checkAndForce(v_u, v_w, G, queue, length_thresh);
                }
            }
        }

        if (!to_visit.empty()) {
            NodeID pick = *to_visit.begin();
            to_visit.erase(pick);
            bool result = explore(G, pick, queue, faces_in_queue, avg_edge_length, length_thresh);
        }
    }
}

/**
 * @brief Build minimum spanning tree (MST)
 *
 * @param mst: [OUT] constructed MST
 * @param p: connection information of the mst
 * @param isEuclidean: if to use Euclidean distance
 * @param vertices: coordinates of the point cloud
 * @param normals: normal of the point cloud
 *
 * @return None
 */
void build_mst(SimpGraph& g, NodeID root,
    RSGraph& out_mst, std::vector<Vector>& normals, std::vector<Point>& vertices) {
    RSGraph temp;
    minimum_spanning_tree(g, root, temp, normals, vertices);

    // Fix strong ambiguous points
    if (!isEuclidean) {
        for (int i = 0; i < temp.m_edges.size(); i++) {
            NodeID source = temp.m_edges[i].source;
            NodeID target = temp.m_edges[i].target;
            Vector normal1 = temp.m_vertices[source].normal;
            normal1 /= normal1.length();
            Vector normal2 = temp.m_vertices[target].normal;
            normal2 /= normal2.length();
            Point pos1 = temp.m_vertices[source].coords;
            Point pos2 = temp.m_vertices[target].coords;
            if (temp.valence(source) >= 2 && temp.valence(target) >= 2)
                continue;
            Vector edge = pos2 - pos1;

            Vector normal_sum = normal1 + normal2;
            float cos_angle = std::abs(dot(edge, normal_sum / normal_sum.length() / edge.length()));
            if (cos_angle > std::cos(10. / 180. * M_PI)) {
                NodeID leaf, parent;
                if (temp.valence(source) == 1) {
                    temp.m_vertices[source].normal = temp.m_vertices[target].normal;
                    parent = target;
                    leaf = source;
                }
                else {
                    temp.m_vertices[target].normal = temp.m_vertices[source].normal;
                    parent = source;
                    leaf = target;
                }

                auto neighbors = temp.neighbors(parent);
                for (auto neighbor : neighbors) {
                    if (temp.m_vertices[neighbor].normal_rep == -1) {
                        temp.m_vertices[neighbor].normal_rep = parent;
                    }
                    else {
                        // Collision!
                        temp.m_vertices[neighbor].normal_rep = -2;
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
        Vector edge = out_mst.m_vertices[temp.m_edges[i].source].coords -
            out_mst.m_vertices[temp.m_edges[i].target].coords;
        float Euclidean_dist = edge.length();
        float projection_dist = cal_proj_dist(edge,
            out_mst.m_vertices[temp.m_edges[i].source].normal,
            out_mst.m_vertices[temp.m_edges[i].target].normal);
        if (std::isnan(projection_dist) || std::isnan(Euclidean_dist))
            std::cout << "debug" << std::endl;

        if(isEuclidean)
            out_mst.add_edge(temp.m_edges[i].source,
                temp.m_edges[i].target, Euclidean_dist);
        else
            out_mst.add_edge(temp.m_edges[i].source,
                temp.m_edges[i].target, projection_dist);
    }
    return;
}

void reset_static() {
    isGTNormal = true;
    isEuclidean = true;
    isFaceLoop = true;
    isDebug = false;
    isNoiseExperiment = false;
    k = 30;
    r = 20.;
    theta = 60.;
    step_thresh = 50;
    exp_genus = -1;
    model_path = "";
    root_path = "";
    model_name = "";
    mode = "";
    recon_timer = RsR_Timer();
    bettiNum_1 = 0;
}
}

namespace HMesh
{

using namespace detail;

bool is_valid_point(const CGLA::Vec3d& p)
{
    return std::isfinite(p[0]) &&
        std::isfinite(p[1]) &&
        std::isfinite(p[2]);
}

bool is_reasonable(const CGLA::Vec3d& p)
{
    const double MAX_VAL = 1e6;

    return std::abs(p[0]) < MAX_VAL &&
        std::abs(p[1]) < MAX_VAL &&
        std::abs(p[2]) < MAX_VAL;
}

void reconstruct_single(HMesh::Manifold& output, std::vector<Vec3d>& org_vertices,
    std::vector<Vec3d>& org_normals, bool use_Euclid_dist, int genus,
    int num_neighbors, int max_neighbor_dist, int max_normal_ang, int max_handle_dist) {
    isEuclidean = use_Euclid_dist;
    exp_genus = genus;
    k = num_neighbors;
    r = max_neighbor_dist;
    theta = max_normal_ang;
    step_thresh = max_handle_dist;

    recon_timer.create("Whole process");
    recon_timer.create("Initialization");
    recon_timer.create("Import pc");
    recon_timer.create("Estimate normals");
    recon_timer.create("Build MST");
    recon_timer.create("Build Rotation System");
    recon_timer.create("algorithm");

    recon_timer.start("Initialization");
    // Estimate normals & orientation & weighted smoothing
    recon_timer.start("Estimate normals");
    std::vector<Point> in_smoothed_v;
    {
        std::vector<Vec3d> clean, clean_normal;
        
        int ii = 0;
        int invalid_num = 0;
        for (auto& p : org_vertices)
        {
            if (!is_valid_point(p) || !is_reasonable(p)) {
                invalid_num++;
                continue;
            }
            clean.push_back(p);
            if (org_normals.size() == org_vertices.size())
                clean_normal.push_back(org_normals[ii]);
            ii++;
        }

        org_vertices = clean;
        org_normals = clean_normal;
        clean_normal.clear();
        clean.clear();
        std::cout << invalid_num << " invalid vertices are removed" << std::endl;

        std::vector<NodeID> indices(org_vertices.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree;
        build_KDTree(kdTree, org_vertices, indices);

        //remove_duplicate_vertices(org_vertices, org_normals, tree_before_remove);

        /*indices.clear();
        indices = std::vector<NodeID>(org_vertices.size());
        std::iota(indices.begin(), indices.end(), 0);
        build_KDTree(kdTree, org_vertices, indices);*/

        float diagonal_length;

        if (isGTNormal) {
            if (!examine_valid_normal(org_normals, org_vertices.size())) {
                isGTNormal = false;
                std::cout << "No input normal / input normal has issue, will estimate normal by myself" << std::endl;
            }
        }

        std::vector<NodeID> zero_normal_id;
        if (!isGTNormal) {
            std::cout << "Estimating normals for a point cloud of size " << org_vertices.size() << std::endl;
            estimate_normal(org_vertices, kdTree, org_normals, zero_normal_id, diagonal_length);
        }

        // Fix zero normal
        if (zero_normal_id.size() > 0) {
            throw std::runtime_error("Zero normal exists!");
        }

        std::cout << "Start first round smoothing ..." << std::endl;
        if (!isEuclidean) {
            weighted_smooth(org_vertices, in_smoothed_v, org_normals, kdTree, diagonal_length);
        }
        else {
            in_smoothed_v = org_vertices;
        }

        Tree temp_tree1;
        build_KDTree(temp_tree1, in_smoothed_v, indices);

        if (!isGTNormal)
            estimate_normal(in_smoothed_v, temp_tree1, org_normals,
                zero_normal_id, diagonal_length);

        // Another round of smoothing
        if (!isEuclidean) {
            std::cout << "Start second round smoothing ..." << std::endl;
            std::vector<Point> temp(in_smoothed_v.begin(), in_smoothed_v.end());
            in_smoothed_v.clear();
            weighted_smooth(temp, in_smoothed_v, org_normals, temp_tree1, diagonal_length);

            if (!isGTNormal) {
                Tree temp_tree2;
                build_KDTree(temp_tree2, in_smoothed_v, indices);

                estimate_normal(in_smoothed_v, temp_tree2, org_normals,
                    zero_normal_id, diagonal_length);
            }
        }
    }
    recon_timer.end("Estimate normals");
    recon_timer.end("Initialization");

    recon_timer.start("algorithm");
    // Find components
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<std::vector<Vector>> component_normals;
    {
        // Build kdTree - CGLA
        std::vector<NodeID> indices(in_smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);

        // Insert number_of_data_points in the tree
        Tree kdTree;
        build_KDTree(kdTree, in_smoothed_v, indices);

        // Correct orientation

        if (!isGTNormal) {
            std::cout << "correct normal orientation" << std::endl;
            correct_normal_orientation(in_smoothed_v, kdTree, org_normals);
        }

        std::cout << "find components" << std::endl;

        find_components(org_vertices, component_vertices, in_smoothed_v,
            component_smoothed_v, org_normals, component_normals, kdTree,
            theta, r);

        in_smoothed_v.clear();
    }

    for (int component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ..." << std::endl;

        isFaceLoop = true;

        std::vector<std::vector<NodeID>> faces;
        std::vector<Point> vertices = component_vertices[component_id];
        std::vector<Vector> normals = component_normals[component_id];
        std::vector<Point> smoothed_v = component_smoothed_v[component_id];

        std::vector<NodeID> indices(smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);

        // Insert number_of_data_points in the tree
        Tree kdTree;
        build_KDTree(kdTree, smoothed_v, indices);

        recon_timer.start("Build MST");

        std::cout << "Init mst" << std::endl;

        // Initial Structure
        RSGraph mst;
        std::vector<m_Edge> full_edges;
        std::vector<m_Edge_length> edge_length;
        std::vector<float> connection_max_length(vertices.size(), 0.);
        std::vector<float> pre_max_length(vertices.size(), 0.);
        mst.isEuclidean = isEuclidean;
        mst.exp_genus = exp_genus;
        SimpGraph g;
        {
            init_graph(smoothed_v, smoothed_v, normals,
                kdTree, g, connection_max_length,
                pre_max_length, theta);

            // Generate MST
            build_mst(g, 0, mst, normals, smoothed_v);

            // Edge arrays and sort
            for (NodeID node : g.node_ids()) {
                for (NodeID node_neighbor : g.neighbors(node)) {
                    if (node < node_neighbor) {
                        Vector edge = smoothed_v[node] - smoothed_v[node_neighbor];
                        double len = edge.length();

                        if (!isEuclidean) {
                            len = cal_proj_dist(edge, normals[node], normals[node_neighbor]);
                        }

                        if (len > pre_max_length[node] ||
                            len > pre_max_length[node_neighbor])
                            continue;
                        edge_length.push_back(m_Edge_length(len, full_edges.size()));
                        full_edges.push_back(m_Edge(node, node_neighbor));
                    }
                }
            }
            std::sort(edge_length.begin(), edge_length.end(), edge_comparator);
        }
        recon_timer.end("Build MST");

        // Initialize face loop label
        mst.etf.reserve(6 * vertices.size() - 11);
        init_face_loop_label(mst);

        // Betti number changes
        std::vector<int> betti_1;

        // Vanilla MST imp
        bettiNum_1 = 0;
        betti_1.push_back(bettiNum_1);
        {
            // Edge connection
            for (int i = 0; i < edge_length.size(); i++) {
                unsigned int edge_idx = edge_length[i].second;
                m_Edge this_edge = full_edges[edge_idx];

                if (mst.find_edge(this_edge.first, this_edge.second) != AMGraph::InvalidEdgeID)
                    continue;

                bool isValid = Vanilla_check(mst, this_edge, kdTree);

                if (isValid) {
                    bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces, kdTree, edge_length[i].first);
                    if (isAdded)
                        betti_1.push_back(bettiNum_1);
                }
            }
        }
        // Create handles & Triangulation
        if (exp_genus != 0) {
            mst.isFinal = true;
            std::vector<NodeID> connected_handle_root;
            connect_handle(smoothed_v, kdTree, mst, connected_handle_root, betti_1);
            isFaceLoop = false;
            triangulate(faces, mst, g, kdTree, isFaceLoop, isEuclidean, connection_max_length, connected_handle_root, betti_1);
        }

        betti_1.push_back(bettiNum_1);
        HMesh::Manifold res;
        // Extract vertex position
        std::vector<double> pos;
        for (int i = 0; i < mst.no_nodes(); i++) {
            pos.push_back(vertices[i][0]);
            pos.push_back(vertices[i][1]);
            pos.push_back(vertices[i][2]);
        }

        std::vector<int> mesh_faces(faces.size(), 3);
        std::vector<int> flattened_face;
        for (int i = 0; i < faces.size(); i++) {
            flattened_face.push_back(faces[i][0]);
            flattened_face.push_back(faces[i][1]);
            flattened_face.push_back(faces[i][2]);
        }

        HMesh::build(res, mst.no_nodes(), &pos[0], faces.size(),
            &mesh_faces[0], &flattened_face[0]);
        output.merge(res);
    }

    recon_timer.end("algorithm");
    std::string line(40, '=');
    std::cout << line << std::endl << std::endl;
    recon_timer.show();
    reset_static();

    return;
}
}