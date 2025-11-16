//
//  comb_quad.cpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 05/03/2023.
//  Copyright © 2023 J. Andreas Bærentzen. All rights reserved.
//

#include "quad_valencify.h"

#include <map>
#include <queue>
#include <algorithm>
#include <random>

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/HMesh.h>

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/HMesh/face_loop.h>

using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Geometry;

struct PQE {
    double priority;
    HalfEdgeID h;
    int time_stamp;
};

bool operator<(const PQE& e0, const PQE& e1) {
    return e0.priority < e1.priority;
}

double ang_weight(const Vec3d& a, const Vec3d& b, const Vec3d& c) {
    double d = dot(a-b, c-b) / (length(a-b)*length(c-b));
    return acos(d);
}

double priority(const HMesh::Manifold& m, HalfEdgeID h) {
    Walker w = m.walker(h);
    
    VertexID v0 = w.vertex();
    Walker wo = w.opp();
    VertexID v1 = wo.vertex();
    
    double l = length(m.pos(v0)-m.pos(v1));
    double ang0 = ang_weight(m.pos(v1), m.pos(v0), m.pos(w.next().vertex())) +
        ang_weight(m.pos(v1), m.pos(v0), m.pos(w.opp().prev().opp().vertex()));
    double ang1 = ang_weight(m.pos(v0), m.pos(v1), m.pos(wo.next().vertex())) +
        ang_weight(m.pos(v0), m.pos(v1), m.pos(wo.opp().prev().opp().vertex()));

    double v0v = valency(m, v0)-4;
    double v1v = valency(m, v1)-4;
    
    if(v0v <= 0 && v1v <= 0)
        return 0;
    if(v0v <= -1 || v1v <= -1)
        return 0;

    return (v0v+v1v)/sqr(ang0+ang1);
}


bool connect_val3(HMesh::Manifold& m) {
    AMGraph3D g;
    Util::AttribVec<AMGraph::EdgeID, HalfEdgeID> gedge2hedge;
    Util::AttribVec<AMGraph::NodeID, VertexID> val3vertex;
    Util::AttribVec<AMGraph::NodeID, FaceID> node2face;
    FaceAttributeVector<AMGraph::NodeID> face2node;

    // Create graph where the faces are the nodes.
    for(auto f: m.faces()) {
        auto n = g.add_node(centre(m, f));
        face2node[f] = n;
        node2face[n] = f;
        val3vertex[n] = InvalidVertexID;
    }
    
    // For each face that contains a valency 3 vertex, store the vertex
    // in a list of valency 3 vertices in the face. If there is one such 
    // vertex for the face, store it in val3vertex for the node.
    vector<pair<VertexID, FaceID>> initial_vertex_face_pairs;
    VertexAttributeVector<int> usage(m.allocated_vertices(), 0);
    for(auto f: m.faces()) {
        VertexID val3_in_face = InvalidVertexID;
        for (auto v: m.incident_vertices(f)) {
            if(valency(m,v)==3)
                if(val3_in_face == InvalidVertexID || usage[v] < usage[val3_in_face]) {
                    val3_in_face = v;
                    usage[v]++;
                }
        }
        if(val3_in_face != InvalidVertexID) {
            AMGraph3D::NodeID n = face2node[f];
            val3vertex[n] = val3_in_face;
            initial_vertex_face_pairs.push_back(make_pair(val3_in_face, f));
        }
    }

    // Connect nodes (faces) if the edge that separates them does not 
    // have a valence 3 vertex at both ends.    
    for(auto h: m.halfedges()) {
        auto w = m.walker(h);
        if(h < w.opp().halfedge())
            if (valency(m, w.vertex())>3 && valency(m, w.opp().vertex()) >3) {
                FaceID f1 = w.face();
                FaceID f2 = w.opp().face();
                auto e = g.connect_nodes(face2node[f1], face2node[f2]);
                gedge2hedge[e] = h;
            }
    }


    for (auto [v_first, f_first]: initial_vertex_face_pairs) {
        BreadthFirstSearch bfs(g);
        bfs.add_init_node(face2node[f_first]);
        while(bfs.Dijkstra_step()) {
            AMGraph::NodeID n_last = bfs.get_last();
            VertexID v_last = val3vertex[n_last];
            if (v_last != InvalidVertexID && v_last != v_first) {
                VertexID v = v_last;
                AMGraph::NodeID n = n_last;
                AMGraph::NodeID np = bfs.pred[n];
                for(; np != AMGraph::InvalidNodeID; np= bfs.pred[n]) {
                    HalfEdgeID h = gedge2hedge[g.find_edge(n, np)];
                    VertexID vs = m.split_edge(h);
                    m.split_face_by_edge(node2face[n], v, vs);
                    v = vs;
                    n = np;
                }
                FaceID f_new =  m.split_face_by_edge(node2face[n], v, v_first);
                if (f_new != InvalidFaceID)
                    return true;
                return false;
            }
        }
    }
    return false;
}

bool connect_same_face_val3(HMesh::Manifold& m) {
    bool did_work = false;
    for(auto f: m.faces()) {
        VertexID v0 = InvalidVertexID;
        VertexID v1 = InvalidVertexID;
        int i0 = 0;
        int i=0;
        int N = no_edges(m, f);
        for(auto v: m.incident_vertices(f)) {
            if (valency(m, v) == 3) {
                if(v0 == InvalidVertexID) {
                    v0 = v;
                    i0 = i;
                }
                else if (i-i0>1 && (i-i0+1) < N){
                    v1 = v;
                    break;
                }
            }
            ++i;
        }
        if (v0 != InvalidVertexID && v1 != InvalidVertexID) {
            if(m.split_face_by_edge(f, v0, v1) != InvalidFaceID)
                did_work = true;
        }
    }
    return did_work;
}

void reduce_valency(HMesh::Manifold& m) {
    
    HalfEdgeAttributeVector<int> time_stamp;
    priority_queue<PQE> Q;
    
    for(auto h: m.halfedges()) {
        auto w = m.walker(h);
        if(h<w.opp().halfedge()) {
            PQE e;
            e.priority = priority(m, h);
            e.time_stamp = time_stamp[h] = 0;
            e.h = h;
            if (e.priority>0)
                Q.push(e);
        }
    }
    
    while(!Q.empty()) {
        auto e = Q.top();
        Q.pop();
        
        if(e.time_stamp == time_stamp[e.h]) {
            Walker w = m.walker(e.h);
            VertexID v0 = w.vertex();
            VertexID v1 = w.opp().vertex();
            m.merge_faces(w.face(), e.h);
            for (auto v : {v0, v1})
                for (auto h : m.incident_halfedges(v)) {
                    auto hmin = m.walker(h).hmin();
                    time_stamp[hmin] = time_stamp[hmin]+1;
                    PQE e;
                    e.priority = priority(m, hmin);
                    if (e.priority>0) {
                        e.time_stamp = time_stamp[hmin];
                        e.h = hmin;
                        Q.push(e);
                    }
                }
        }
    }
}

void quad_valencify(HMesh::Manifold& m) {
    reduce_valency(m);
    while (connect_same_face_val3(m));
    while (connect_val3(m));
    for(auto v: m.vertices())
        if (valency(m, v) == 3)
            cout << "VAL 3 LEFT!!!" << endl;
}

struct VertexPath {
    vector<VertexID> vertices;
    vector<HalfEdgeID> edges;
    int priority = 0;
    
    void join(const VertexPath& other) {
        vertices.insert(end(vertices), begin(other.vertices), end(other.vertices));
        edges.insert(end(edges), begin(other.edges), end(other.edges));
    }
    
    bool simple() {
        std::sort(vertices.begin(), vertices.end());
        return std::adjacent_find(vertices.begin(), vertices.end()) == vertices.end();
    }
    
    void compute_priority(const vector<VertexID>& orig_verts) {
        for (auto v: orig_verts) {
            if(find(begin(vertices), end(vertices), v) != end(vertices))
                priority += 1;
        }
    }
    
    size_t size() const { return vertices.size();}
    
    bool operator<(const VertexPath& other) const {
        return priority < other.priority;
    }
    
    void print() {
        cout << "Path: " << vertices.size() << " vertices " << edges.size() << " edges " << " pri: " << priority << endl;
    }
};

vector<VertexPath> trace_vertex_paths(const Manifold& m,
                                      const vector<VertexID>& o_verts) {
    
    auto trace = [&](const DijkstraOutput& d_out, VertexID v_start) {
        VertexPath v_path;
        for(auto v = v_start; d_out.pred[v] != InvalidVertexID; v = d_out.pred[v]) {
            v_path.vertices.push_back(v);
            v_path.edges.push_back(m.walker(d_out.pred_edge[v]).hmin());
        }
        return v_path;
    };
    
    vector<VertexPath> v_paths;
    for (auto v0: o_verts) {
        DijkstraOutput d_out = Dijkstra(m, v0);
        
        for (auto v_leaf: d_out.leaves) {
            for (auto h: m.incident_halfedges(v_leaf)) {
                auto w = m.walker(h);
                auto vn = w.vertex();
                if (vn != d_out.pred[v_leaf]) {
                    auto v_path = trace(d_out, v_leaf);
                    v_path.join(trace(d_out, vn));
                    v_path.edges.push_back(w.hmin());
                    v_path.vertices.push_back(v0);
                    if(v_path.size() >= 6 && v_path.simple()) {
                        v_path.compute_priority(o_verts);
                        v_paths.push_back(v_path);
                    }
                }
            }
        }
    }
    return v_paths;
}


void quad_valencify_no_crossings(HMesh::Manifold& m,
                                 HMesh::VertexAttributeVector<CGLA::Vec3f>& vcol,
                                 HMesh::HalfEdgeAttributeVector<CGLA::Vec3f>& hcol) {
    
    vector<VertexID> orig_verts;
    for (auto v: m.vertices())
        orig_verts.push_back(v);
    
    vector<HalfEdgeID> orig_edges;
    for (auto h: m.halfedges())
        if (h == m.walker(h).hmin())
            orig_edges.push_back(h);
    
    for (auto f:m.faces())
        m.split_face_by_vertex(f);
    
    for (auto h: orig_edges) {
        Walker w = m.walker(h);
        FaceID f0 = w.face();
        VertexID v0 = w.next().vertex();
        FaceID f1 = w.opp().face();
        VertexID v1 = w.opp().next().vertex();
        VertexID v = m.split_edge(h);
        m.split_face_by_edge(f0, v, v0);
        m.split_face_by_edge(f1, v, v1);
    }
    
    auto paths = trace_vertex_paths(m, orig_verts);
    VertexAttributeVector<int> v_touched(m.no_vertices(), 0);
    HalfEdgeAttributeVector<int> h_touched(m.no_halfedges(),0);
    priority_queue<VertexPath> Q(begin(paths), end(paths));
    while (! Q.empty()) {
        auto path = Q.top();
        Q.pop();
        bool valid = true;
        for(auto v: path.vertices) {
            if(v_touched[v]==2) {
                valid = false;
                break;
            }
        }
        if (valid) for(auto h: path.edges) {
            if(h_touched[h] > 0) {
                valid = false;
                break;
            }
        }
        
        bool all_done = false;
        if (valid) {
            all_done = true;
            for (auto v: path.vertices)
                v_touched[v] += 1;
            for (auto h: path.edges)
                h_touched[h] += 1;
            
            
            for (auto v: orig_verts) {
                if (v_touched[v] < 2) {
                    all_done = false;
                    break;
                }
            }
        }
        if (all_done) {
            cout << "All done!!" << endl;
            break;
        }
        
//        if(valid) {
//            cout << "adding just one " << endl;
//            break;
//        }
    }
    if(Q.empty())
        cout << "Drained Q" << endl;
    
    for (auto v: m.vertices()) {
        cout << "Vertex: " << v.index << " touched " << v_touched[v];
        if (find(begin(orig_verts), end(orig_verts), v) != orig_verts.end())
            cout << " is orig ";
        cout << endl;
        vcol[v] = Vec3f(v_touched[v]==2,v_touched[v]==1,v_touched[v]==0);
    }
    
    for (auto h: m.halfedges()) {
        auto w = m.walker(h);
        if (h_touched[w.hmin()] > 1)
            cout << "Edge: " << w.hmin().index << " touched: " << h_touched[w.hmin()] << endl;
        hcol[h] = Vec3f(h_touched[w.hmin()]==1,0,0);
    }
    
    return;
    
    for (auto v: m.vertices()) {
        if(m.in_use(v) && v_touched[v]==0) {
            if (find(begin(orig_verts), end(orig_verts), v) != orig_verts.end())
                cout << "bad bad news ... this is an original vertex: " << v_touched[v] << endl;
            cout << "merging one ring of vertex: " << v.index << endl;
            m.merge_one_ring(v);
//            m.remove_vertex(v);
//            close_holes(m);
        }
    }
    for(auto _h: m.halfedges())
        if (m.in_use(_h)) {
            auto w = m.walker(_h);
            auto h = w.hmin();
            if(h_touched[h]==0) {
                cout << "Merging: " << h.index;
                bool success = m.merge_faces(w.face(), w.halfedge());
                cout << " success: " << success << endl;
            }
    }
}

// Function to hash an integer to a color
Vec3f hashToColor(int value) {
    // Prime numbers for hashing
    const int prime1 = 73856093;
    const int prime2 = 19349663;
    const int prime3 = 83492791;
    
    value += 421312;

    // Generate color components using prime numbers
    float r = (value * prime1) % 256 / 255.0f;
    float g = (value * prime2) % 256 / 255.0f;
    float b = (value * prime3) % 256 / 255.0f;
    return Vec3f (r, g, b);
}

void index_preserving_cc_split(HMesh::Manifold& m) {
    VertexAttributeVector<int> touched(m.no_vertices(), -1);
    // Mark all original vertices
    for (auto v: m.vertices())
        touched[v] = 1;
    
    // Split all of the original edges by inserting
    // a vertex on them.
    vector<HalfEdgeID> orig_edges;
    for (auto h: m.halfedges())
        if (h == m.walker(h).hmin())
            orig_edges.push_back(h);
    for (auto h: orig_edges)
        m.split_edge(h);
    
    // Split all vertices by inserting a vertex in the middle
    vector fvec(begin(m.faces()), end(m.faces()));
    for(auto f: fvec) {
        auto vnew = m.split_face_by_vertex(f);
        touched[vnew] = 2;
    }

    // Remove edges from the face points to original vertices.
    for(auto h: m.halfedges())
        if (m.in_use(h)) {
            auto w = m.walker(h);
            if(touched[w.vertex()] == 1 && touched[w.opp().vertex()] == 2)
                m.merge_faces(w.face(), h);
        }
    
}

void quad_valencify_cc(HMesh::Manifold& m_orig,
                       HMesh::VertexAttributeVector<CGLA::Vec3f>& vcol,
                       HMesh::HalfEdgeAttributeVector<CGLA::Vec3f>& hcol,
                       HMesh::FaceAttributeVector<CGLA::Vec3f>& fcol) {
    bool all_v4 = true;
    for(auto v: m_orig.vertices()) {
        if (valency(m_orig, v) != 4)
            all_v4 = false;
    }
    
    if (all_v4)
        return;
    
    
    Manifold m = m_orig;
    VertexSet orig_verts = m.vertices();
   
    index_preserving_cc_split(m);
    
    FaceAttributeVector<VertexID> face_to_orig_vert(m.no_faces(),InvalidVertexID);
    map<VertexID,int> face_count;
    for(auto v: m.vertices()) {
        if (orig_verts.contains(v)) {
            face_count[v] = 0;
            for (auto f: m.incident_faces(v)) {
                face_to_orig_vert[f] = v;
                fcol[f] = hashToColor(v.index);
                face_count[v] += 1;
            }
        }
    }
//    vector<HalfEdgeID> hvec(begin(m.halfedges()), end(m.halfedges()));
//    std::random_device rd;
//    std::mt19937 gen {rd()};
//    shuffle(begin(hvec), end(hvec), gen);
    bool did_work = true;
    while (did_work) {
        did_work = false;
        auto loops = find_face_loops(m);
        
        face_loops_compute_contained_area(m, loops);
        sort(begin(loops), end(loops), [](const FaceLoop& f1, const FaceLoop& f2) {
            return f1.interior_area<f2.interior_area;
        });
        
        for (auto& fl: loops) {
            map<VertexID,int> loop_face_count;
            for (auto v: orig_verts)
                loop_face_count[v] = 0;
            bool removable = true;

            for (auto h_fl: fl.hvec) {
                Walker w = m.walker(h_fl);
                loop_face_count[face_to_orig_vert[w.face()]] += 1;
                if (valency(m, w.vertex()) + valency(m, w.opp().vertex()) == 6) {
                    removable = false;
                    break;
                }
            }

            if (removable)
                for(auto v: orig_verts)
                    if  (loop_face_count[v]==face_count[v])
                        removable = false;
            
            if (removable) {
                remove_face_loop(m, fl);
                for(auto v: orig_verts)
                    face_count[v]-=loop_face_count[v];
                did_work = true;
                break;
            }
        }
    }
    
    for (auto v: m.vertices())
        m.pos(v).normalize();
    
    Manifold m_dual;
    for (auto v: m.vertices()) {
        // add_face cannot be called with owned views
        auto incident = m.incident_faces(v);
        auto pos = std::views::all(incident)
        | std::views::transform([&m](FaceID fn) { return centre(m, fn); });
        m_dual.add_face(pos);
    }
    stitch_mesh(m_dual,1e-10);
    for (auto v: m_dual.vertices())
        m_dual.pos(v).normalize();

    m_dual.cleanup();
    
    m_orig = m_dual;
    
}
