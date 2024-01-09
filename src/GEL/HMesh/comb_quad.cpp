//
//  comb_quad.cpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 05/03/2023.
//  Copyright © 2023 J. Andreas Bærentzen. All rights reserved.
//

#include "comb_quad.h"
#include <array>
#include <queue>

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/Range.h>

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_util.h>

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

    return  l*(v0v+v1v)/ sqr(max(ang0,ang1));
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
}


