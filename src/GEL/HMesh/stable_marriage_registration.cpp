//
//  stable_marriage_matching.cpp
//  MeshEditE
//
//  Created by Jakob Andreas Bærentzen on 16/01/2025.
//  Copyright © 2025 J. Andreas Bærentzen. All rights reserved.
//

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include "stable_marriage_matching.h"
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/build_bbtree.h>

using namespace HMesh;
using namespace CGLA;
using namespace Geometry;
//
//double stable_marriage_match(HMesh::Manifold& m, const HMesh::Manifold& m_ref) {
//    // The inputs m and m_ref are polygonal meshes. They are stored using the GEL api:
//    // https://github.com/janba/GEL
//    
//    // Make a point vector, P, with all vertices of m. Generate points on m
//    // and add them to P until |P| == m_ref.no_vertices().
//    // For each point, p, of P, store associated normal and closest vertex of m.
//    // The closest vertex of m is chosen from among those incident on p.
//    
//    // For each point in P, compute pseudo distance to all vertices of m_ref.
//    // Pseudo distance should be computed as d * (1+epsilon-dot(n, n_ref))
//    // In other words if two points are close (d is small) or the dot product
//    // between the normal n and the normal n_ref from the reference mesh is
//    // close to 1, the pseudo distance is small. For each vertex of m_ref make
//    // and sort increasing order a vector of pseudo distances and indices of points in P.
//    
//    
//    // Now, we loop, computing a stable marriage of vertices in m_ref and points in P.
//    // For each unassigned vertex, v, in m_ref:
//    //   if v does not have an associated point, loop over P in order of increasing pseudo distance until
//    //   you find a point, p, in P for which the
//    //   vertex. If there is no currently assigned vertex, v is simply assigne. Otherwise if the pseudo distance from v
//    //   to p is smaller than the pseudo distance from v' to p where v' is the currently assigned vertex, then v' is
//    //   unassigned and v is assigned to p.
//    // This loop is repeated while unassigned vertices remain. It will terminate since the number of vertices in m_ref
//    // is equal to the number of points in P.
//    // This is the Gale–Shapley algorithm.
//    
//    // For each point, p, of P, compute the vector to the assigned vertex, v, of m_ref. for every vertex of
//    // m average these vectors for all associated points in P. Move vertex by adding this averaged vector
//    // to its position.
//    // return sum of pseudo distances between points in P and their final assigned vertices.
//    
//    return 0;
//}
//
//
//#include "HMesh/Manifold.h"
//#include "HMesh/VertexAttributeVector.h"
//
//using namespace HMesh;

using namespace std;

Vec3d rand_point_in_face(const Manifold& m, FaceID f) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double rand_sum = 0;
    Vec3d p(0);
    for (auto v: m.incident_vertices(f)) {
        double random_value = dis(gen);
        rand_sum += random_value;
        p += m.pos(v) * random_value;
    }
    return p / rand_sum;
}

struct P_type {
    // Create a point vector P with all vertices of m
    Vec3d pos;
    Vec3d norm;
    FaceID face;
    VertexID v_ref;
};

VertexAttributeVector<set<int>> trace_to_surface(const Manifold& m_ref, const vector<P_type>& P) {
    
    cout << "Building AABB tree" << endl;
    AABBTree T;
    build_AABBTree(m_ref, T);
    int N = P.size();
    double ael = average_edge_length(m_ref);
    vector<double> dists(N,0);
    vector<Vec3d> pts(N);
    cout << "Tracing to surface" << endl;
    for (int i=0; i<N; ++i) {
        Vec3f p(P[i].pos);
        Vec3f dir(P[i].norm);
        auto d = T.compute_signed_distance(p);
        if (abs(d)>ael) {
            p -= d * dir;
            dists[i] += abs(d);
        }
        pts[i] = Vec3d(p);
    }
    cout << "Building KD tree of reference vertices" << endl;

    KDTree<Vec3d, VertexID> vertex_tree;
    for (auto v: m_ref.vertices())
        vertex_tree.insert(m_ref.pos(v), v);
    vertex_tree.build();
    
    cout << "creating index sets for each vertex" << endl;
    VertexAttributeVector<set<int>> vertex_pt_sets;
    for (int i=0;i<N;++i) {
        KDTree<Vec3d, VertexID>::ScalarType d = 10.0 * ael;
        Vec3d k;
        VertexID v0;
        if (vertex_tree.closest_point(pts[i], d, k, v0)) {
            vertex_pt_sets[v0].insert(i);
            queue<VertexID> vertex_queue;
            vertex_queue.push(v0);
            while(not vertex_queue.empty()) {
                VertexID v = vertex_queue.front();
                vertex_queue.pop();
                for (auto w: m_ref.incident_vertices(v)) {
                    if (not vertex_pt_sets[w].contains(i) and length(m_ref.pos(v0)-m_ref.pos(w))<max(dists[i], 2*ael) ) {
                        vertex_queue.push(w);
                        vertex_pt_sets[w].insert(i);
                    }
                }
            }
        }
        cout << endl;
    }
    cout << "bye from trace" << endl;
    return vertex_pt_sets;
}

double stable_marriage_match(Manifold& m, Manifold& m_ref, int N) {
    
    m.cleanup();
    m_ref.cleanup();
    
    // Define epsilon for pseudo distance calculation
    double search_dist;
    float r;
    Vec3d c;
    bsphere(m_ref, c, r);
    search_dist = DBL_MAX;

    
    // Generate points on m until |P| == m_ref.no_vertices()
    FaceAttributeVector<double> face_areas;
    double total_area = 0.0;
    for (auto f : m.faces()) {
        double a = area(m,f);
        face_areas[f] = a;
        total_area += a;
    }

    std::vector<double> cumulative_areas;
    double cumulative_sum = 0.0;
    for (auto f : m.faces()) {
        cumulative_sum += face_areas[f] / total_area;
        cumulative_areas.push_back(cumulative_sum);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);


    std::vector<P_type> P;
    FaceAttributeVector<vector<int>> face2P_idx;
    while (P.size() < m_ref.no_vertices()) {
        double random_value = dis(gen);
        auto it = std::lower_bound(cumulative_areas.begin(), cumulative_areas.end(), random_value);
        int face_index = std::distance(cumulative_areas.begin(), it);
        FaceID f = FaceID(face_index);
        Vec3d p = rand_point_in_face(m, f);
        face2P_idx[f].push_back(P.size());
        P.push_back({p, normal(m,f), f, InvalidVertexID});
    }
    cout << "Computing distances ... " << endl;
    // For each point in P, compute pseudo distance to all vertices of m_ref
    
    auto vertex_pt_sets = trace_to_surface(m_ref, P);
    VertexAttributeVector<std::vector<std::pair<double, int>>> pseudo_distances;
    for (auto v_ref : m_ref.vertices()) {
        const Vec3d& p_ref = m_ref.pos(v_ref);
        for (int i: vertex_pt_sets[v_ref]) {
            double d = length(P[i].pos-p_ref)*(1.1-dot(P[i].norm, normal(m_ref, v_ref)));
            pseudo_distances[v_ref].emplace_back(d, i);
        }
    }
    cout << "on to sorting " << endl;

    // Sort pseudo distances for each vertex in m_ref
    int zero_cnt = 0;
    for (auto v_ref: m_ref.vertices()) {
        if (pseudo_distances[v_ref].size() > 0)
            std::sort(pseudo_distances[v_ref].begin(), pseudo_distances[v_ref].end());
        else zero_cnt += 1;
    }
    
    cout << zero_cnt << " vertices lacked matches?!" << endl;

    cout << "Stable marriage matching ..." << endl;
    // Stable marriage algorithm
    VertexAttributeVector<int> assignment(m_ref.no_vertices(), -1);
    std::vector<bool> assigned(P.size(), false);
    bool unassigned_vertices = true;
    bool did_work = true;
    while (unassigned_vertices and did_work) {
        unassigned_vertices = false;
        did_work = false;
        for (auto v_ref: m_ref.vertices()) {
            if (assignment[v_ref] == -1) {
                unassigned_vertices = true;
                for (const auto& [pseudo_distance, p_idx] : pseudo_distances[v_ref]) {
                    if (!assigned[p_idx]) {
                        assignment[v_ref] = p_idx;
                        assigned[p_idx] = true;
                        P[p_idx].v_ref = v_ref;
                        did_work = true;
                        break;
                    } else {
                        VertexID current_v_ref = P[p_idx].v_ref;
                        if (pseudo_distance < pseudo_distances[current_v_ref][0].first) {
                            assignment[current_v_ref] = -1;
                            assignment[v_ref] = p_idx;
                            P[p_idx].v_ref = v_ref;
                            did_work = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    cout << "Move vertices ..." << endl;

    for (auto v : m.vertices()) {
        Vec3d avg_vector(0, 0, 0);
        int count = 0;
        
        for (auto f: m.incident_faces(v)) {
            for (int p_idx: face2P_idx[f]) {
                const auto& pt = P[p_idx];
                if (pt.v_ref != InvalidVertexID) {
                    avg_vector += m_ref.pos(pt.v_ref) - pt.pos;
                    count += 1;
                }
            }
        }
        if (count > 0) {
            avg_vector /= count;
            m.pos(v) += avg_vector;
        }
        else {
            cout << "vertex " << v.index << " came up empty" << endl;
        }
    }

    return 0;//total_pseudo_distance;
}
