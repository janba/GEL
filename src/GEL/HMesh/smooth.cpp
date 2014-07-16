/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <thread>
#include "smooth.h"

#include <future>
#include <vector>
#include <algorithm>
#include "../CGLA/Mat3x3d.h"
#include "../CGLA/Vec3d.h"
#include "../CGLA/Quatd.h"
#include "../Util/Timer.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    
    
    class range {
    public:
        class iterator {
            friend class range;
        public:
            long int operator *() const { return i_; }
            const iterator &operator ++() { ++i_; return *this; }
            iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }
            
            bool operator ==(const iterator &other) const { return i_ == other.i_; }
            bool operator !=(const iterator &other) const { return i_ != other.i_; }
            
        protected:
            iterator(long int start) : i_ (start) { }
            
        private:
            unsigned long i_;
        };
        
        iterator begin() const { return begin_; }
        iterator end() const { return end_; }
        range(long int  begin, long int end) : begin_(begin), end_(end) {}
    private:
        iterator begin_;
        iterator end_;
    };
    
    void for_each_vertex(Manifold& m, function<void(VertexID)> f) { for(auto v : m.vertices()) f(v); }
    void for_each_face(Manifold& m, function<void(FaceID)> f) { for(auto p : m.faces()) f(p); }
    void for_each_halfedge(Manifold& m, function<void(HalfEdgeID)> f) { for(auto h : m.halfedges()) f(h); }

    
    inline Vec3d laplacian(const Manifold& m, VertexID v)
    {
        Vec3d p(0);
        int n = circulate_vertex_ccw(m, v, (std::function<void(VertexID)>)[&](VertexID v){ p += m.pos(v); });
        return p / n - m.pos(v);
    }
    
    
    int CORES = 8;
    
    void laplacian_smooth_example(Manifold& m)
    {
        VertexAttributeVector<Vec3d> new_pos = m.positions_attribute_vector();
        for(VertexID v : m.vertices())
            if(!boundary(m, v))
            {
                Vec3d L(0);
                int n = circulate_vertex_ccw(m, v, (std::function<void(VertexID)>)[&](VertexID v){ L += m.pos(v); });
                new_pos[v] = L/n;
            }
        m.positions_attribute_vector() = new_pos;
    
    }
    
    void laplacian_smooth0(Manifold& m,float weight, int max_iter)
    {
        for(int iter=0;iter<max_iter; ++iter) {
            VertexAttributeVector<Vec3d> L_attr(m.no_vertices());
            
            for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
                if(!boundary(m, *v))
                    L_attr[*v] =laplacian(m, *v);
            for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                if(!boundary(m, *v))
                    m.pos(*v) += weight*L_attr[*v];
            }
        }
    }
    
    void laplacian_smooth1(Manifold& m,float weight, int max_iter)
    {
        for(int iter=0;iter<max_iter; ++iter) {
            VertexAttributeVector<Vec3d> L_attr(m.no_vertices());
            
            for(VertexID v : m.vertices())
                if(!boundary(m, v))
                    L_attr[v] =laplacian(m, v);
            for(VertexID v : m.vertices()){
                if(!boundary(m, v))
                    m.pos(v) += weight*L_attr[v];
            }
        }
    }
    
    void laplacian_smooth2(Manifold& m,float weight, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter=0;iter<max_iter; ++iter) {
            for(auto v : m.vertices())
                if(!boundary(m, v))
                    new_pos[v] = weight*laplacian(m, v)+m.pos(v);
            m.positions_attribute_vector() = new_pos;
        }
    }


    
    void laplacian_smooth3(Manifold& m, float weight, int max_iter)
    {
        for(int iter=0;iter<max_iter; ++iter) {
            auto new_pos = m.positions_attribute_vector();
            for(auto v : m.vertices())
                if(!boundary(m, v))
                    new_pos[v] = weight*laplacian(m, v)+m.pos(v);
            m.positions_attribute_vector() = move(new_pos);
        }
    }

    void laplacian_smooth4(Manifold& m,float weight, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter=0;iter<max_iter; ++iter) {
            for(VertexID v : m.vertices())
                if(!boundary(m, v))
                    new_pos[v] = weight*laplacian(m, v)+m.pos(v);
            swap(m.positions_attribute_vector(),new_pos);
        }
    }
    
    void laplacian_smooth4_5(Manifold& m,float weight, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter=0;iter<max_iter; ++iter) {
            for_each_vertex(m, [&](VertexID v) {new_pos[v] = weight*laplacian(m, v)+m.pos(v);});
            swap(m.positions_attribute_vector(),new_pos);
        }
    }


    
    void laplacian_smooth5(Manifold& m, float weight, int max_iter)
    {
        for(int iter=0;iter<max_iter; ++iter) {
            auto new_pos = m.positions_attribute_vector();
            vector<thread> t_vec;
            for(auto v : m.vertices())
                if(!boundary(m, v))
                    t_vec.push_back(thread([&](VertexID vid){
                        if(!boundary(m, vid))
                            new_pos[vid] = weight*laplacian(m, vid)+ m.pos(vid);},v));
            for(int i=0;i<t_vec.size();++i)
                t_vec[i].join();
            m.positions_attribute_vector() = move(new_pos);
        }
    }
    
    
    
    inline void laplacian_smooth_vertex(Manifold& m, const vector<VertexID>& vids,
                                        VertexAttributeVector<Vec3d>& new_pos,
                                        float weight){
        for(VertexID v: vids)
            new_pos[v] = m.pos(v)+weight*laplacian(m, v);
    }
    
    
    void laplacian_smooth6(Manifold& m, float weight, int max_iter)
    {
        vector<vector<VertexID>> vertex_ids(CORES);
        auto batch_size = m.no_vertices()/CORES;
        int cnt = 0;
        for_each_vertex(m, [&](VertexID v) {
            if (!boundary(m, v))
                vertex_ids[(cnt++/batch_size)%CORES].push_back(v);
        });
        vector<thread> t_vec(CORES);
        VertexAttributeVector<Vec3d> new_pos = m.positions_attribute_vector();
        for(int iter=0;iter<max_iter; ++iter) {
            for(int thread_no=0;thread_no<CORES;++thread_no)
                t_vec[thread_no] = thread(laplacian_smooth_vertex,
                                          ref(m), ref(vertex_ids[thread_no]),
                                          ref(new_pos), weight);
            for(int thread_no=0;thread_no<CORES;++thread_no)
                t_vec[thread_no].join();
            swap(m.positions_attribute_vector(), new_pos);
        }
        
    }
    
    
    typedef std::vector<std::vector<VertexID>> VertexIDBatches;

//    template<typename  T>
//    void for_each_vertex_parallel(int no_threads, const VertexIDBatches& batches, const T& f) {
//        vector<future<void>> f_vec(no_threads);
//        for(auto t : range(0, no_threads))
//            f_vec[t] = async(launch::async, f, ref(batches[t]));
//    }
    
    
    template<typename  T>
    void for_each_vertex_parallel(int no_threads, const VertexIDBatches& batches, const T& f) {
        vector<thread> t_vec(no_threads);
        for(auto t : range(0, no_threads))
            t_vec[t] = thread(f, ref(batches[t]));
        for(auto t : range(0, no_threads))
            t_vec[t].join();
    }
    
    VertexIDBatches batch_vertices(Manifold& m) {
        VertexIDBatches vertex_ids(CORES);
        auto batch_size = m.no_vertices()/CORES;
        int cnt = 0;
        for_each_vertex(m, [&](VertexID v) {
            //if (!boundary(m, v))
                vertex_ids[(cnt++/batch_size)%CORES].push_back(v);
        });
        return vertex_ids;
    }

    void laplacian_smooth7(Manifold& m, float weight, int max_iter)
    {
        auto vertex_ids = batch_vertices(m);
        auto new_pos = m.positions_attribute_vector();
        auto f = [&](const vector<VertexID>& vids) {
            for(VertexID v: vids)
                new_pos[v] = m.pos(v)+weight*laplacian(m, v);
        };

        for(auto _ : range(0, max_iter)) {
            for_each_vertex_parallel(CORES, vertex_ids, f);
            swap(m.positions_attribute_vector(), new_pos);
        }
    }

 


    void laplacian_smooth(Manifold& m,float weight, int max_iter)
    {
//        laplacian_smooth_example(m);
        laplacian_smooth7(m,weight,max_iter);

//        Util::Timer tim;
//        cout << "Method 0: ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth0(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        cout << "Method 1: ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth1(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        cout << "Method 2: ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth2(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        cout << "Method 3: ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth3(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        cout << "Method 4: ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth4(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        cout << "Method 5: ";
//        for(int i=0;i<3;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth5(m2, 1, 1000);
//            cout << " sec: "<< tim.get_secs() << ",";
//        }
//        cout << endl;
//
//        CORES = 2;
//        cout << "Method 6 (2): ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth6(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//        
//        CORES = 4;
//        cout << "Method 7 (4): ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth6(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;
//
//        CORES = 8;
//        cout << "Method 7 (8): ";
//        for(int i=0;i<5;++i) {
//            Manifold m2 = m;
//            tim.start();
//            laplacian_smooth6(m2, 1, 1000);
//            cout << tim.get_secs() << ", ";
//        }
//        cout << endl;


    }

    

    
    
    void taubin_smooth(Manifold& m, int max_iter)
    {
        auto new_pos = m.positions_attribute_vector();
        for(int iter = 0; iter < 2*max_iter; ++iter) {
            for(VertexID v : m.vertices())
                if(!boundary(m, v))
                    new_pos[v] = (iter%2 == 0 ? 0.5 : -0.52) * laplacian(m, v) + m.pos(v);
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
            
            FaceAttributeVector<Vec3d> filtered_norms(m.allocated_faces());
            
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
                    cout << "iters " << sub_iter << endl;
                    break;
                }
            }
        }
    }
    
    void TAL_smoothing(Manifold& m, float w, int max_iter)
    {
        for(int iter=0;iter<max_iter;++iter) {
            VertexAttributeVector<float> vertex_areas;
            VertexAttributeVector<Vec3d> laplacians;
            
            for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            {
                vertex_areas[*vid] = 0;
                for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
                    if(w.face() != InvalidFaceID)
                        vertex_areas[*vid] += area(m, w.face());
            }
            
            for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            {
                laplacians[*vid] = Vec3d(0);
                double weight_sum = 0.0;
                if(boundary(m, *vid))
                {
                    double angle_sum = 0;
                    for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
                    {
                        if (w.face() != InvalidFaceID)
                        {
                            Vec3d vec_a = normalize(m.pos(w.vertex()) - m.pos(*vid));
                            Vec3d vec_b = normalize(m.pos(w.circulate_vertex_ccw().vertex()) -
                                                    m.pos(*vid));
                            angle_sum += acos(max(-1.0,min(1.0,dot(vec_a,vec_b))));
                        }
                        if(boundary(m,w.vertex()))
                        {
                            laplacians[*vid] += m.pos(w.vertex()) - m.pos(*vid);
                            weight_sum += 1.0;
                        }
                    }
                    laplacians[*vid] /= weight_sum;
                    laplacians[*vid] *= exp(-3.0*sqr(max(0.0, M_PI-angle_sum)));
                }
                else
                {
                    for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
                    {
                        float weight = vertex_areas[w.vertex()];
                        Vec3d l = m.pos(w.vertex()) - m.pos(*vid);
                        laplacians[*vid] +=  weight * l;
                        weight_sum += weight;
                    }
                    laplacians[*vid] /= weight_sum;
                    //                Vec3d n = normal(m, *vid);
                    //                if(sqr_length(n)>0.9)
                    //                    laplacians[*vid] -= n * dot(n, laplacians[*vid]);
                }
            }
            for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
                m.pos(*vid) += w*laplacians[*vid];
        }
    }

    
    
}
