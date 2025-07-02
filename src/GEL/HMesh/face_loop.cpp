//
//  gem.cpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 24/11/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#include "face_loop.h"

#include <array>
#include <queue>

//#include "LogMap.h"
#include <GEL/CGLA/CGLA.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/Range.h>

using namespace Util;
using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace GLGraphics;
using namespace HMesh;

void create_hex(Manifold& m, vector<Vec3d>& pts, double scale) {
    Vec3d avg_p(0);
    for(const auto& p: pts)
        avg_p += p;
    avg_p /= pts.size();
    
    for(auto& p: pts)
        p = scale*(p-avg_p) + avg_p;
    
    m.add_face({pts[0],pts[1],pts[2],pts[3]});
    m.add_face({pts[7],pts[6],pts[5],pts[4]});
    m.add_face({pts[4],pts[5],pts[1],pts[0]});
    m.add_face({pts[5],pts[6],pts[2],pts[1]});
    m.add_face({pts[6],pts[7],pts[3],pts[2]});
    m.add_face({pts[7],pts[4],pts[0],pts[3]});
}

double geodesic_curvature(Manifold& m, Walker& w) {
    Vec3d v = m.pos(w.vertex()) - m.pos(w.opp().vertex());
    Walker wn = w.next().next();
    Vec3d vn = m.pos(wn.opp().vertex())-m.pos(wn.vertex());
    return acos(min(1.0,max(-1.0,dot(v,vn)/(length(v)*length(vn)))));
};


double aspect(Manifold& m, FaceID f) {
    double lmin=DBL_MAX;
    double lmax=1e-30;
    for (auto h: m.incident_halfedges(f)) {
        double l = length(m, h);
        lmin = min(l, lmin);
        lmax = max(l, lmax);
    }
    return lmin/lmax;
}



FaceLoop trace_face_loop(Manifold& m, HalfEdgeAttributeVector<int>& touched, HalfEdgeID h) {
    auto comp_valency_imbalance = [](Manifold& m, Walker& w) {
        int val1 = valency(m, w.vertex());
        int val2 = valency(m, w.opp().vertex());
        
        if(val1 < 4 && val2 > 4)
            return val2-val1;
        if(val2 < 4 && val1 > 4)
            return val1-val2;
        if(val1 < 4 && val2 < 4)
            return -1;
        return 0;
    };
    
    Walker w = m.walker(h);
    
    FaceLoop loop;
    FaceAttributeVector<int> path_faces_touched(0);
    Vec3d avg_edge(0);
    Vec3d center(0);
    double avg_len = 0;
    do {
        touched[w.halfedge()] = 1;
        loop.area += area(m, w.face());
        if(path_faces_touched[w.face()] == 1)
            loop.double_cross_faces.insert(w.face());
        path_faces_touched[w.face()] = 1;
        loop.hvec.push_back(w.halfedge());
        Vec3d p0 = m.pos(w.vertex());
        Vec3d p1 = m.pos(w.opp().vertex());
        Vec3d edg = (p0 - p1);
        double edg_len = edg.length();
        center += p0 + p1;
        avg_edge += edg;
        avg_len += edg_len;
        int imbalance = comp_valency_imbalance(m, w);
        loop.valency_imbalance += imbalance;
        loop.integral_geodesic_curvature += geodesic_curvature(m,w);
        loop.min_len = min(loop.min_len, edg_len);
        double a = aspect(m, w.face());
        loop.min_aspect = min(loop.min_aspect, a);
        loop.mean_aspect += a;
        w = w.next().next().opp();
    }
    while(w.halfedge() != h && !touched[w.halfedge()]);
    loop.cylindricity = length(avg_edge)/avg_len;
    loop.avg_edge = avg_edge;
    loop.avg_len = avg_len / loop.hvec.size();
    loop.center = center / (2.0 * loop.hvec.size());
    loop.integral_geodesic_curvature /= loop.area;
    loop.valency_imbalance /= loop.area;
    loop.mean_aspect /= loop.hvec.size();
//    loop.cylindricity;// /= loop.area;
    
    return loop;
}

vector<FaceLoop> find_face_loops(Manifold& m)
{
    vector<FaceLoop> face_loops;
    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
    int loop_no = 0;
    for(auto h: m.halfedges())
    {
        if(!touched[h]){
            face_loops.push_back(trace_face_loop(m, touched, h));
            FaceLoop& loop = face_loops.back();
            loop.id = ++loop_no;
        }
    }
    return face_loops;
}

void face_loops_compute_contained_area(Manifold& m, vector<FaceLoop>& face_loops)
{
    FaceAttributeVector<int> face_status(m.no_faces(),0);
    for(FaceLoop& l: face_loops) {
        if(l.double_cross_faces.empty()) {
            
            // Mark all faces in the loop with loop id.
            for(auto h: l.hvec)
                face_status[m.walker(h).face()] = l.id;

            queue<FaceID> fq;
            for(auto h: l.hvec) {
                Walker w = m.walker(h);
                FaceID f = w.prev().opp().face();
                if(face_status[f]!= l.id)
                    fq.push(f);
            }
            
            while(!fq.empty())
            {
                auto f = fq.front();
                fq.pop();
                if(l.interior.count(f)==0)
                {
                    l.interior_faces += 1;
                    if (f==InvalidFaceID)
                        cout << "f was invalid" << endl;
                    l.interior_area += sqrt(area(m, f));
                    l.interior.insert(f);
                    circulate_face_ccw(m, f, [&](FaceID fn){
                        if(l.interior.count(fn)==0 && face_status[fn] != l.id)
                            fq.push(fn);
                    });
                }
                
            }
        }
        else l.interior_area = DBL_MAX;
    }
}

void collapse_face_loops(Manifold& m, vector<FaceLoop>& face_loops)
{
    FaceAttributeVector<int> face_status(m.no_faces(),0);
    VertexAttributeVector<int> touched(m.no_vertices(),0);
    for(FaceLoop& l: face_loops) {
        if(l.double_cross_faces.empty()) {
            for(auto h: l.hvec)
                face_status[m.walker(h).face()] = l.id;
            
            queue<FaceID> fq;
            for(auto h: l.hvec) {
                Walker w = m.walker(h);
                FaceID f = w.prev().opp().face();
                if(face_status[f]!= l.id)
                    fq.push(f);
            }
            
            VertexSet vs;
            for(FaceID f: l.interior)
                circulate_face_ccw(m, f, [&](VertexID v){
                    if(!touched[v]) {
                        vs.insert(v);
                        touched[v] = 1;
                    }
                });
            
            Vec3d barycenter(0);
            int n=0;
            for(auto v: vs)
            {
                barycenter += m.pos(v);
                ++n;
            }
            barycenter /= n;
            for(auto v: vs)
                m.pos(v) = 0.9*barycenter+0.1*m.pos(v);
        }
    }
}

int collapse_double_crossed_quads(HMesh::Manifold& m)
{
    auto face_loops = find_face_loops(m);
    double DE_min = DBL_MAX;
    VertexID va_min;
    VertexID vb_min;
    FaceID f_min = InvalidFaceID;
    for(auto& fl: face_loops ) {
        if(!fl.double_cross_faces.empty()) {
            for(auto f: fl.double_cross_faces)
                if(m.in_use(f))
                {
                    Walker w = m.walker(f);
                    VertexID v0 = w.vertex();
                    VertexID v1 = w.next().vertex();
                    VertexID v2 = w.next().next().vertex();
                    VertexID v3 = w.next().next().next().vertex();
                    
                    vector<int> val = {valency(m, v0),valency(m, v1), valency(m, v2), valency(m, v3)};
                    double E_bef = 0;
                    for (auto v : val)
                        E_bef += sqr(v-4.0);
                    
                    double E_after_02 = sqr(val[0]+val[2]-2-4) + sqr(val[1]-1-4) + sqr(val[3]-1-4);
                    double E_after_13 = sqr(val[1]+val[3]-2-4) + sqr(val[0]-1-4) + sqr(val[2]-1-4);
                    
                    double DE02 = E_after_02 - E_bef;
                    double DE13 = E_after_13 - E_bef;
                    double DE;
                    VertexID va,vb;
                    if(DE02 < DE13) {
                        va = v0;
                        vb = v2;
                        DE = DE02;
                    }
                    else {
                        va = v1;
                        vb = v3;
                        DE = DE13;
                    }
                    if (DE < DE_min) {
                        DE_min = DE;
                        va_min = va;
                        vb_min = vb;
                        f_min = f;
                    }
                }
        }
    }
    if(f_min != InvalidFaceID) {
        FaceID f2 = m.split_face_by_edge(f_min, va_min, vb_min);
        if (f2 != InvalidFaceID) {
            Walker w = m.walker(f2);
            if(precond_collapse_edge(m, w.halfedge())) {
                m.collapse_edge(w.halfedge(),true);
                return 1;
            }
        }
    }
    return 0;
}


int split_double_crossed_quads(HMesh::Manifold& m)
{
    int cnt=0;
    auto face_loops = find_face_loops(m);
    FaceSet faces_to_split;
    for(auto& fl: face_loops )
        faces_to_split.insert(begin(fl.double_cross_faces), end(fl.double_cross_faces));
                
    for(auto f: faces_to_split)
    {
        Walker w = m.walker(f); VertexID v0 = w.vertex();
        VertexID v1 = w.next().vertex();
        VertexID v2 = w.next().next().vertex();
        VertexID v3 = w.next().next().next().vertex();
        if(sqr_length(m.pos(v0)-m.pos(v2)) > sqr_length(m.pos(v1)-m.pos(v3))) {
            swap(v0,v1);
            swap(v2,v3);
        }
        auto fnew = m.split_face_by_edge(f, v0, v2);
        m.split_edge(m.walker(fnew).halfedge());
        cnt += 1;
    }
    return cnt;
}





FaceSet extrude_along_edge_loop(Manifold& m, const HalfEdgeSet& hset)
{
    vector<pair<HalfEdgeID, HalfEdgeID>> h_pairs;
    
    for(auto h: hset) {
        auto w = m.walker(h);
        h_pairs.push_back(make_pair(w.halfedge(), w.opp().halfedge()));
    }
    
    // Slit edges
    for(auto h_in: hset) {
        Walker w_in = m.walker(h_in);
        VertexID v = w_in.vertex();
        for(auto h_out: hset) {
            Walker w_out = m.walker(h_out);
            if(w_out.opp().vertex() == v) {
                auto vnew = m.slit_vertex(v, h_in, h_out);
                break;
            }
        }
    }
    
    VertexAttributeVector<int> cluster_id(m.allocated_vertices(),-1);
    int CLUSTER_CNT = 0;
    auto assign_cluster_id = [&](VertexID v) {
        if(cluster_id[v]==-1)
            cluster_id[v]= ++CLUSTER_CNT;
    };
    
    FaceSet new_faces;
    // Make new faces
    for(auto h_pair: h_pairs)
    {
        Walker w1 = m.walker(h_pair.first);
        Walker w2 = m.walker(h_pair.second);
        assign_cluster_id(w1.vertex());
        assign_cluster_id(w1.opp().vertex());
        assign_cluster_id(w2.vertex());
        assign_cluster_id(w2.opp().vertex());

        auto&& pts = {
            m.pos(w1.vertex()),
            m.pos(w1.opp().vertex()),
            m.pos(w2.vertex()),
            m.pos(w2.opp().vertex())
        };
        
        FaceID f = m.add_face(pts);
        new_faces.insert(f);
        Walker w =  m.walker(f);
        
        cluster_id[w.vertex()] = cluster_id[w1.vertex()];
        w = w.next();
        cluster_id[w.vertex()] = cluster_id[w1.opp().vertex()];
        w = w.next();
        cluster_id[w.vertex()] = cluster_id[w2.vertex()];
        w = w.next();
        cluster_id[w.vertex()] = cluster_id[w2.opp().vertex()];
    }
    
    // Stitch
    stitch_mesh(m, cluster_id);
    
    return new_faces;
}


FaceSet extrude_face_set(Manifold& m, const FaceSet& face_set)
{
    HalfEdgeSet hset;
    
    for(auto f: face_set) {
        circulate_face_ccw(m, f, [&](Walker& w){
            if(face_set.find(w.opp().face()) == face_set.end()) {
                hset.insert(w.halfedge());
            }
        });
    }
    return extrude_along_edge_loop(m, hset);
}

FaceSet extrude_halfedge_set(Manifold& m, HalfEdgeSet& halfedge_set)
{
    vector<pair<HalfEdgeID, HalfEdgeID>> h_pairs;
    HalfEdgeSet hset;
    
    for(auto h: halfedge_set)
    {
        Walker w = m.walker(h);
        hset.insert(w.halfedge());
        hset.insert(w.opp().halfedge());
        h_pairs.push_back(make_pair(w.halfedge(), w.opp().halfedge()));
    }
    halfedge_set.clear();
    
    int CLUSTER_CNT = 0;
    VertexAttributeVector<int> orig_id;
    for(auto v: m.vertices())
        orig_id[v] = ++CLUSTER_CNT;
    
    VertexSet new_verts;
    // Slit edges
    for(auto h_in: hset)
        for(auto h_out: hset) {
            Walker w_in = m.walker(h_in);
            Walker w_out = m.walker(h_out);
            VertexID v = w_in.vertex();
            new_verts.insert(w_in.opp().vertex());
            new_verts.insert(v);
            new_verts.insert(w_out.vertex());
            if(w_out.opp().vertex() == v &&
               w_in.halfedge() != w_out.opp().halfedge())
            {
                VertexID v_new = m.slit_vertex(v, h_in, h_out);
                if(v_new != InvalidVertexID) {
                    orig_id[v_new] = orig_id[v];
                    break;
                }
            }
        }
    
    
    VertexAttributeVector<int> cluster_id(m.allocated_vertices(),-1);
    auto assign_cluster_id = [&](VertexID v) {
        if(cluster_id[v]==-1)
            cluster_id[v]= ++CLUSTER_CNT;
    };
    // Make new faces
    FaceSet new_faces;
    for(auto h_pair: h_pairs)
    {
        Walker w1 = m.walker(h_pair.first);
        Walker w2 = m.walker(h_pair.second);
        assign_cluster_id(w1.vertex());
        assign_cluster_id(w1.opp().vertex());
        assign_cluster_id(w2.vertex());
        assign_cluster_id(w2.opp().vertex());
        
        auto&& ptsa = {
            m.pos(w1.vertex()),
            m.pos(w1.opp().vertex()),
            m.pos(w1.opp().vertex()),
            m.pos(w1.vertex())
        };
        
        FaceID fa = m.add_face(ptsa);
        Walker wa = m.walker(fa);
        new_faces.insert(fa);

        cluster_id[wa.vertex()] = cluster_id[w1.vertex()];
        wa = wa.next();
        cluster_id[wa.vertex()] = cluster_id[w1.opp().vertex()];
        wa = wa.next();
        cluster_id[wa.vertex()] = orig_id[w1.opp().vertex()];
        new_verts.insert(wa.vertex());
        wa = wa.next();
        halfedge_set.insert(wa.halfedge());
        new_verts.insert(wa.vertex());
        cluster_id[wa.vertex()] = orig_id[w1.vertex()];
        
        auto&& ptsb =
        {
            m.pos(w2.opp().vertex()),
            m.pos(w2.vertex()),
            m.pos(w2.vertex()),
            m.pos(w2.opp().vertex())
        };
        
        FaceID fb = m.add_face(ptsb);
        Walker wb = m.walker(fb);
        new_faces.insert(fb);

        cluster_id[wb.vertex()] = orig_id[w2.opp().vertex()];
        new_verts.insert(wb.vertex());
        wb = wb.next();
        cluster_id[wb.vertex()] = orig_id[w2.vertex()];
        new_verts.insert(wb.vertex());
        wb = wb.next();
        cluster_id[wb.vertex()] = cluster_id[w2.vertex()];
        wb = wb.next();
        cluster_id[wb.vertex()] = cluster_id[w2.opp().vertex()];
        
    }
    
    // Stitch
    stitch_mesh(m, cluster_id);
    
    for (int iter=0;iter<1;++iter) {
        VertexAttributeVector<Vec3d> new_pos;
        for(auto v: new_verts)
            if(m.in_use(v))
            {
                double cnt = 0.0;
                new_pos[v] = Vec3d(0);
                circulate_vertex_ccw(m, v,[&](VertexID vn){
                    double w = 1.0;
                    new_pos[v] += w*m.pos(vn);
                    cnt += w;
                });
                new_pos[v] /= cnt;
            }
        for(auto v: new_verts)
            if(m.in_use(v))
                m.pos(v) = new_pos[v];
    }
    
    return new_faces;
}


vector<FaceID> create_box(Manifold& m, const Vec3d& pos, const Mat3x3d& _R, double sz)
{
    double h = 0.5;
    const auto face_points = std::array<Vec3d const, 4> {Vec3d(0,-h,-h),Vec3d(0,h,-h),Vec3d(0,h,h),Vec3d(0,-h,h)};
    
    Mat3x3d R = _R;
    double det = determinant(R);
    if(abs(det))
        if(det<0) {
            Mat3x3d M = identity_Mat3x3d();
            M[2][2] = -1;
            R = R * M;
        }
    
    vector<FaceID> fvec;
    for(int axis : {0,1,2})
        for(int sgn : {-1,1}) {
            // maybe better to just use a simple array here?
            auto pts = std::views::all(face_points) |
                std::views::transform([&](Vec3d _p) {
                    Vec3d p(0);
                    p[(0+axis)%3] += _p[0] + h * sgn;
                    p[(1+axis)%3] += sgn * _p[1];
                    p[(2+axis)%3] += _p[2];
                    return (R*(p*sz))+pos;
                });
            fvec.push_back(m.add_face(pts));
        }
    return fvec;
}


bool remove_face_loop(Manifold& m, const FaceLoop& l, bool average) {
    
    VertexAttributeVector<int> cluster_id(m.allocated_vertices(), -1);
    vector<VertexSet> cluster;
    int cid_N = 0;
    
    // The block of code below makes a list of faces to remove but it also
    // clusters vertices such that we can stitch the mesh when the faces have
    // been removed.
    vector<FaceID> fvec;
    for(auto h: l.hvec) {
        Walker w = m.walker(h);
        VertexID v0 = w.vertex();
        VertexID v1 = w.opp().vertex();
        assert(m.in_use(v0));
        assert(m.in_use(v1));
        assert(m.in_use(h));

        int cid_max = max(cluster_id[v0],cluster_id[v1]);
        if(cid_max == -1) {
            cid_max = cid_N;
            cid_N += 1;
            VertexSet c;
            cluster_id[v0] = cid_max;
            cluster_id[v1] = cid_max;
            c.insert(v0);
            c.insert(v1);
            cluster.push_back(c);
        }
        else {
            if(cluster_id[v0]==-1) {
                cluster_id[v0] = cid_max;
                cluster[cid_max].insert(v0);
            }
            else if(cluster_id[v1]==-1) {
                cluster_id[v1] = cid_max;
                cluster[cid_max].insert(v1);
            } else {
                auto cid_min = min(cluster_id[v0],cluster_id[v1]);
                for (auto v : cluster[cid_max])
                    cluster_id[v] = cid_min;
                cluster[cid_min].insert(cluster[cid_max].begin(),cluster[cid_max].end());
            }
        }
        fvec.push_back(w.face());
    }
    
    // Set the position of vertices on the edge of the face loop.
    // We normally see this is an inverse extrusion, so the vertex on the edge
    // that is contracted is assigned the position of the vertex on the side towards
    // which we contract. However, in some cases, multiple vertices are involved.
    // In these cases, we contract to the average position.
    for(auto h: l.hvec) {
        Walker w = m.walker(h);
        VertexID v0 = w.opp().vertex();
        int cid = cluster_id[v0];
        if(!average && cluster[cid].size() == 2)
            m.pos(v0) = m.pos(w.vertex());
        else {
            Vec3d p(0);
            for(auto v : cluster[cid])
                p += m.pos(v);
            p /= cluster[cid].size();
            m.pos(v0) = p;
            m.pos(w.vertex()) = p;
        }
    }
    
    //Remove faces.
    for(auto f: fvec)
        if(m.in_use(f))
            m.remove_face(f);
    
    // Stitch up mesh
    int unstitched = stitch_mesh(m, cluster_id);
    
    // The conditional below is entered only if there is a bug.
    // We should always be able to stitch.
    if (unstitched) {
        cout << "some unstitched: " << unstitched << endl;
        return false;
    }
    return true;
}

void kill_face_loop(Manifold& m) {
    
    auto loops = find_face_loops(m);
    
    face_loops_compute_contained_area(m, loops);
    sort(begin(loops), end(loops), [](const FaceLoop& f1, const FaceLoop& f2) {
        return f1.interior_area<f2.interior_area;
    });
    
    // Visit all face loops
    auto& l  = loops.front();
    cout << "Interior area " << l.interior_area << endl;
    auto& hvec_min = l.hvec;
    // Find the boundary vertices and then the interior vertices
    // the interior vertices are used below when we flatten the interior.
    VertexSet boundary_vertices;
    for(auto h:hvec_min)
    {
        Walker w = m.walker(h);
        boundary_vertices.insert(w.opp().vertex());
        boundary_vertices.insert(w.vertex());
    }
    VertexSet interior_vertices;
    for(FaceID f: l.interior)
        circulate_face_ccw(m, f, [&](VertexID v){
            if(!boundary_vertices.count(v))
                interior_vertices.insert(v);
        });
    
    // Remove the actual loop...
    remove_face_loop(m, l);

    // Finally smooth the patch inside region bounded by the killed loop.
    for(int iter=0;iter<10;++iter)
    {
        for(auto v: interior_vertices)
            if(m.in_use(v))
            {
                Vec3d avg_pos(0);
                avg_pos /= circulate_vertex_ccw(m,v, [&](VertexID vn){avg_pos += m.pos(vn);});
                m.pos(v) = avg_pos;
            }
    }
}

void kill_degenerate_face_loops(Manifold& m, double thresh) {
    double avg_edge_len = average_edge_length(m);
    bool did_work = true;
    int cnt = 0;
    VertexAttributeVector<int> touched(0);
    while (did_work) {
        did_work = false;
        auto loops = find_face_loops(m);
        sort(begin(loops), end(loops), [](const FaceLoop& f1, const FaceLoop& f2) {
            return f1.min_len<f2.min_len;
        });
        if (loops[0].min_len < thresh*avg_edge_len) {
            bool allowable = true;
            for (auto h: loops[0].hvec) {
                auto w = m.walker(h);
                if ( (valency(m, w.vertex())<4 and valency(m, w.opp().vertex())<4) or
                    (touched[w.vertex()] or touched[w.opp().vertex()]) ) {
                    allowable = false;
                    break;
                }
            }
            if (allowable) {
                did_work = true;
                for (auto h: loops[0].hvec) {
                    auto w = m.walker(h);
                    touched[w.vertex()] = 1;
                    touched[w.opp().vertex()] = 1;
                }
                remove_face_loop(m, loops[0], true);
                ++cnt;
            }
        }
    }
//    cout << "Killed " << cnt << " face loops " << endl;
    m.cleanup();
}


int refine_loop(Manifold& m) {
    
    auto loops = find_face_loops(m);
    sort(loops.begin(), loops.end(), [&](FaceLoop& f0, FaceLoop& f1){
        return f0.avg_len > f1.avg_len;
    });
    
    int work = 0;
    for(auto l: loops)
        if(l.double_cross_faces.size()==0)
        {
            HalfEdgeSet hset;
            bool dismiss = false;
            for(auto h: l.hvec) {
                Walker w =m.walker(h).next();
                if(hset.count(w.opp().halfedge()) == 0)
                    hset.insert(w.halfedge());
                else
                    dismiss = true;
            }
            if(dismiss)
                continue;
            extrude_along_edge_loop(m, hset);
            ++work;
            break;
        }
    return work;
}

Geometry::AMGraph3D face_loop_skeleton(HMesh::Manifold& m) {
    
    AMGraph3D graph;
    
    auto loops = find_face_loops(m);
    
    int f_no = m.no_faces();
    sort(loops.begin(), loops.end(), [&](FaceLoop& f0, FaceLoop& f1){
        return f0.cylindricity > f1.cylindricity;
    });
    
    HalfEdgeAttributeVector<AMGraph::NodeID> h2n_map(m.allocated_halfedges(), AMGraph::InvalidNodeID);
    FaceAttributeVector<int> touched(m.allocated_faces(),0);
    for(auto i: Util::Range(0,loops.size()))
    {
        auto& l = loops[i];
        if(l.double_cross_faces.empty()) {
            bool face_loop_intersects = false;
            for (auto h: l.hvec)
            {
                FaceID f = m.walker(h).face();
                if(touched[f])
                    face_loop_intersects = true;
            }
            
            if(!face_loop_intersects) {
                AMGraph::NodeID n0,n1;
                
                // Create graph nodes corresponding to the face loop, we kill.
                n0 = graph.add_node(Vec3d(0));
                n1 = graph.add_node(Vec3d(0));
                
                auto e = graph.connect_nodes(n0, n1);
                graph.edge_color[e] = Vec3f(0,0,0.6);
                
                // Move vertices to new position corresponding to the loop being removed.
                for(auto h: l.hvec) {
                    Walker w = m.walker(h);
                    FaceID f = w.face();
                    touched[f] = 1;
                    DebugRenderer::face_colors[f] = Vec3f(1.0-l.cylindricity,0.3,l.cylindricity);
                    VertexID v0 = w.vertex();
                    VertexID v1 = w.opp().vertex();
                    graph.pos[n0] += m.pos(v0);
                    graph.pos[n1] += m.pos(v1);
                    
                    HalfEdgeID h0 = w.next().halfedge();
                    HalfEdgeID h1 = w.prev().halfedge();
                    h2n_map[h0] = n0;
                    h2n_map[h1] = n1;
                }
                
                // Finalize position of end points of skeletal edge.
                graph.pos[n0] /= l.hvec.size();
                graph.pos[n1] /= l.hvec.size();
            }
        }
    }
    
    for(auto f0: m.faces())
        if(touched[f0] == 0)
        {
            
            AMGraph::NodeID patch_node = graph.add_node(Vec3d(0));
            graph.node_color[patch_node] = Vec3f(1,0,0);
            queue<FaceID> fq;
            fq.push(f0);
            touched[f0] = 2;
            int i=0;
            while(!fq.empty())
            {
                FaceID f = fq.front();
                fq.pop();
                circulate_face_ccw(m, f, [&](Walker w){
                    FaceID fo = w.opp().face();
                    if(touched[fo] == 0) {
                        fq.push(fo);
                        touched[fo] = 2;
                    }
                    else if(touched[fo]==1) {
                        h2n_map[w.halfedge()] = patch_node;
                        AMGraph::NodeID no = h2n_map[w.opp().halfedge()];
                        if(graph.valid_node_id(no)) {
                            graph.pos[patch_node] += graph.pos[no];
                            ++i;
                        }
                        else
                            cout << "Bad node id" << endl;
                    }
                });
            }
            graph.pos[patch_node] /= i;
        }
    
    
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        auto e0 = h2n_map[h];
        auto e1 = h2n_map[w.opp().halfedge()];
        if(e0 != AMGraph::InvalidNodeID && e1 != AMGraph::InvalidNodeID) {
            auto e = graph.connect_nodes(e0, e1);
            if(e != AMGraph::InvalidEdgeID)
                graph.edge_color[e] = Vec3f(0.2,0.9,0.4);
        }
    }
    return graph;
}
