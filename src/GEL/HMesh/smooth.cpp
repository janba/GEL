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
            double a_left  = acos(min(0.99, max(-0.99, d_left)));
            double a_right = acos(min(0.99, max(-0.99, d_right)));
            
//            double w = sin(a_left + a_right) / (0.001+sin(a_left)*sin(a_right));
            double w = cos(a_left)/sin(a_left) + cos(a_right)/sin(a_right);
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

    void taubin_smooth_cot(Manifold& m, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter = 0; iter < 2*max_iter; ++iter) {
            for(VertexID v : m.vertices())
                new_pos[v] = (iter%2 == 0 ? +0.5 : -0.52) * cot_laplacian(m, v) + m.pos(v);
            swap(m.positions_attribute_vector(), new_pos);
        }
    }

vector<FaceID> face_neighbourhood(Manifold& m, FaceID f)
{
    vector<FaceID> nbrs;
    Walker w = m.walker(*m.incident_halfedges(f).begin());
    int N = no_edges(m, f);
    for(int i=0; i<N; ++i) {
        do {
            FaceID fn = w.face();
            if (fn != InvalidFaceID)
                nbrs.push_back(fn);
            w = w.next().opp();
        }
        while (w.face() != f);
        w = w.opp().next().opp();
    }
    return nbrs;
}

    
    void anisotropic_smooth(HMesh::Manifold& m, int max_iter, double sharpness)
    {
        for(int iter = 0;iter<max_iter; ++iter)
        {
            FaceAttributeVector<Vec3d> filtered_norms;
            FaceAttributeVector<Vec3d> norms;
            for (FaceID f: m.faces())
                norms[f] = normal(m, f);
            
            for (FaceID f: m.faces()) {
                Vec3d n0 = norms[f];
                Vec3d n_new(0);
                for (FaceID fn: face_neighbourhood(m, f)) {
                    Vec3d n = norms[fn];
                    double w_a = exp(-sharpness*acos(max(-1.0,min(1.0,dot(n,n0))))/M_PI);
                    n_new += area(m, fn) * w_a * n;
                }
                filtered_norms[f] = normalize(n_new);
            }
            
            VertexAttributeVector<Vec3d> vertex_positions(m.allocated_vertices(), Vec3d(0));
            VertexAttributeVector<int> count(m.allocated_vertices(), 0);
            for (auto f: m.faces()) {
                Vec3d n_o = normal(m, f);
                Quatd q;
                q.make_rot(n_o, filtered_norms[f]);
                Mat3x3d M = q.get_Mat3x3d();
                Vec3d c = barycenter(m, f);
                for (auto v: m.incident_vertices(f)) {
                    vertex_positions[v] += c + M * (m.pos(v)-c);
                    count[v] += 1;
                }
            }
            for(VertexID v : m.vertices())
                m.pos(v) = vertex_positions[v] / double(count[v]);
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

    

     void regularize_quads(HMesh::Manifold& m, float weight, float shrink, int max_iter) {
        for(int iter = 0; iter < max_iter; ++iter) {
            auto new_pos = VertexAttributeVector<Vec3d>(m.allocated_vertices(), Vec3d(0));
            auto cnt = VertexAttributeVector<int>(m.allocated_vertices(), 0);
            for (auto f: m.faces()) {
                vector<Vec3d> corners;
                vector<VertexID> corner_vids;
                for (auto v: m.incident_vertices(f)) {
                    corners.push_back(m.pos(v));
                    corner_vids.push_back(v);
                }

                if(corners.size() != 4)
                    continue;
                for (auto v: corner_vids)
                    cnt[v] += 1;

                Vec3d c = corners[0]+corners[1]+corners[2]+corners[3];
                c /= 4.0;

                Vec3d a = corners[2]-corners[0];
                Vec3d b = corners[3]-corners[1];
                double l = (1.0-shrink)*0.25 * (length(a) + length(b));
                a.normalize();
                b.normalize();
                new_pos[corner_vids[0]] += c - l * a;
                new_pos[corner_vids[1]] += c - l * b;
                new_pos[corner_vids[2]] += c + l * a;
                new_pos[corner_vids[3]] += c + l * b;
            }

            for (auto v: m.vertices()) {
                m.pos(v) = weight * new_pos[v] / cnt[v] + (1-weight)*m.pos(v);  
            }
        }


     }
    
}
