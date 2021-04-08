/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <thread>
#include <GEL/HMesh/smooth.h>

#include <future>
#include <vector>
#include <algorithm>
#include <GEL/CGLA/Mat3x3d.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Quatd.h>
#include <GEL/Util/Timer.h>
#include <GEL/Util/Range.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>

namespace HMesh
{
    using namespace std;
    using namespace CGLA;

    void for_each_vertex(Manifold& m, function<void(VertexID)> f) { for(auto v : m.vertices()) f(v); }
    void for_each_face(Manifold& m, function<void(FaceID)> f) { for(auto p : m.faces()) f(p); }
    void for_each_halfedge(Manifold& m, function<void(HalfEdgeID)> f) { for(auto h : m.halfedges()) f(h); }

    
    int CORES = 8;
    
    typedef std::vector<std::vector<VertexID>> VertexIDBatches;
    
    template<typename  T>
    void for_each_vertex_parallel(int no_threads, const VertexIDBatches& batches, const T& f) {
        vector<thread> t_vec(no_threads);
        for(auto t : Util::Range(0, no_threads))
            t_vec[t] = thread(f, ref(batches[t]));
        for(auto t : Util::Range(0, no_threads))
            t_vec[t].join();
    }
    
    VertexIDBatches batch_vertices(Manifold& m) {
        VertexIDBatches vertex_ids(CORES);
        auto batch_size = m.no_vertices()/CORES;
        int cnt = 0;
        for_each_vertex(m, [&](VertexID v) {
                vertex_ids[(cnt++/batch_size)%CORES].push_back(v);
        });
        return vertex_ids;
    }

    void laplacian_smooth(Manifold& m, float weight, int max_iter)
    {
        auto vertex_ids = batch_vertices(m);
        auto new_pos = m.positions_attribute_vector();
        auto f = [&](const vector<VertexID>& vids) {
            for(VertexID v: vids)
                new_pos[v] = m.pos(v)+weight*laplacian(m, v);
        };
        for(int i=0; i < max_iter; ++i) {
            for_each_vertex_parallel(CORES, vertex_ids, f);
            swap(m.positions_attribute_vector(), new_pos);
        }
    }


    CGLA::Vec3d cot_laplacian(const Manifold& m, VertexID v)
    {
        CGLA::Vec3d p(0);
        Vec3d vertex = m.pos(v);
        double w_sum=0.0;
        circulate_vertex_ccw(m, v, [&](Walker wv){
            Vec3d nbr(m.pos(wv.vertex()));
            Vec3d left(m.pos(wv.next().vertex()));
            Vec3d right(m.pos(wv.opp().prev().opp().vertex()));
            
            double d_left = dot(cond_normalize(nbr-left),
                                cond_normalize(vertex-left));
            double d_right = dot(cond_normalize(nbr-right),
                                 cond_normalize(vertex-right));
            double a_left  = acos(min(1.0, max(-1.0, d_left)));
            double a_right = acos(min(1.0, max(-1.0, d_right)));
            
            double w = sin(a_left + a_right) / (1e-10+sin(a_left)*sin(a_right));
            p += w * nbr;
            w_sum += w;
        });
        if(w_sum<1e-20 || std::isnan(p[0])  || std::isnan(p[1]) || std::isnan(p[2]))
            return Vec3d(0);
        return p / w_sum - m.pos(v);
    }
    
    
    void taubin_smooth(Manifold& m, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter = 0; iter < 2*max_iter; ++iter) {
            for(VertexID v : m.vertices())
                new_pos[v] = (iter%2 == 0 ? +0.5 : -0.52) * laplacian(m, v) + m.pos(v);
            swap(m.positions_attribute_vector(), new_pos);
        }
    }
    
    void face_neighbourhood(Manifold& m, FaceID f, vector<FaceID>& nbrs)
    {
        FaceAttributeVector<int> touched(m.allocated_faces(), -1);
        touched[f] = 1;
        nbrs.push_back(f);
        for(Walker wf = m.walker(f); !wf.full_circle(); wf = wf.circulate_face_cw()){
            for(Walker wv = m.walker(wf.vertex()); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
                FaceID fn = wv.face();
                if(fn != InvalidFaceID && touched[fn] != touched[f]){
                    nbrs.push_back(fn);
                    touched[fn] = 1;
                }
            }
        }
    }
    
    
    
    Vec3d fvm_filtered_normal(Manifold& m, FaceID f)
    {
        const float sigma = .1f;
        
        vector<FaceID> nbrs;
        face_neighbourhood(m, f, nbrs);
        float min_dist_sum=1e32f;
        long int median=-1;
        
        vector<Vec3d> normals(nbrs.size());
        for(size_t i=0;i<nbrs.size();++i)
        {
            normals[i] = normal(m, nbrs[i]);
            float dist_sum = 0;
            for(size_t j=0;j<nbrs.size(); ++j)
                dist_sum += 1.0f - dot(normals[i], normals[j]);
            if(dist_sum < min_dist_sum)
            {
                min_dist_sum = dist_sum;
                median = i;
            }
        }
        assert(median != -1);
        Vec3d median_norm = normals[median];
        Vec3d avg_norm(0);
        for(size_t i=0;i<nbrs.size();++i)
        {
            float w = exp((dot(median_norm, normals[i])-1)/sigma);
            if(w<1e-2) w = 0;
            avg_norm += w*normals[i];
        }
        return normalize(avg_norm);
    }
    Vec3d bilateral_filtered_normal(Manifold& m, FaceID f, double avg_len)
    {
        vector<FaceID> nbrs;
        face_neighbourhood(m, f, nbrs);
        Vec3d p0 = centre(m, f);
        Vec3d n0 = normal(m, f);
        vector<Vec3d> normals(nbrs.size());
        vector<Vec3d> pos(nbrs.size());
        Vec3d fn(0);
        for(FaceID nbr : nbrs)
        {
            Vec3d n = normal(m, nbr);
            Vec3d p = centre(m, nbr);
            double w_a = exp(-acos(max(-1.0,min(1.0,dot(n,n0))))/(M_PI/32.0));
            double w_s = exp(-length(p-p0)/avg_len);
            
            fn += area(m, nbr)* w_a * w_s * n;
            
        }
        return normalize(fn);
    }
    
    void anisotropic_smooth(HMesh::Manifold& m, int max_iter, NormalSmoothMethod nsm)
    {
        double avg_len=0;
        for(HalfEdgeID hid : m.halfedges())
            avg_len += length(m, hid);
        avg_len /= 2.0;
        for(int iter = 0;iter<max_iter; ++iter)
        {
            
            FaceAttributeVector<Vec3d> filtered_norms;
            
            for(FaceID f: m.faces()){
                filtered_norms[f] = (nsm == BILATERAL_NORMAL_SMOOTH)?
                bilateral_filtered_normal(m, f, avg_len):
                fvm_filtered_normal(m, f);
            }
            
            VertexAttributeVector<Vec3d> vertex_positions(m.allocated_vertices(), Vec3d(0));
            VertexAttributeVector<int> count(m.allocated_vertices(), 0);
            for(int sub_iter=0;sub_iter<100;++sub_iter)
            {
                for(HalfEdgeID hid : m.halfedges())
                {
                    Walker w = m.walker(hid);
                    FaceID f = w.face();
                    
                    if(f != InvalidFaceID){
                        VertexID v = w.vertex();
                        Vec3d dir = m.pos(w.opp().vertex()) - m.pos(v);
                        Vec3d n = filtered_norms[f];
                        vertex_positions[v] += m.pos(v) + 0.5 * n * dot(n, dir);
                        count[v] += 1;
                    }
                    
                }
                double max_move= 0;
                for(VertexID v : m.vertices())
                {
                    Vec3d npos = vertex_positions[v] / double(count[v]);
                    double move = sqr_length(npos - m.pos(v));
                    if(move > max_move)
                        max_move = move;
                    
                    m.pos(v) = npos;
                }
                if(max_move<sqr(1e-8*avg_len))
                {
//                    cout << "iters " << sub_iter << endl;
                    break;
                }
            }
        }
    }

    void TAL_smoothing(Manifold& m, float w, int max_iter)
    {
        for(int iter=0;iter<max_iter;++iter) {
            FaceAttributeVector<Vec3d> face_normal;
            FaceAttributeVector<double> face_area;
            for(auto f: m.faces()) {
                face_normal[f] = normal(m, f);
                face_area[f] = area(m, f);
            }
            
            VertexAttributeVector<float> vertex_area;
            VertexAttributeVector<Vec3d> L;
            VertexAttributeVector<Vec3d> norm;
            

            for(auto v: m.vertices())
            {
                vertex_area[v] = 0;
                norm[v] = Vec3d(0);
                for (auto f: m.incident_faces(v))
                    if(m.in_use(f)) {
                        vertex_area[v] += face_area[f];
                        norm[v] += face_normal[f] * face_area[f];
                    }
                norm[v].cond_normalize();
            }
            
            for(auto v: m.vertices())
            {
                L[v] = Vec3d(0);
                if(!boundary(m, v)) {
                    double weight_sum = 0.0;
                    for(auto vn: m.incident_vertices(v))
                    {
                        float weight = vertex_area[vn];
                        Vec3d l = m.pos(vn) - m.pos(v);
                        L[v] +=  weight * l;
                        weight_sum += weight;
                    }
                    L[v] /= weight_sum;
                    Vec3d n = norm[v];
                    if(sqr_length(n)>0.9)
                        L[v] -= n * dot(n, L[v]);
                }
            }
            for(auto v: m.vertices())
                m.pos(v) += w*L[v];
        }
    }

    
    
}
