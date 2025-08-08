/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <queue>
#include <iostream>
#include <thread>
#include <algorithm>

#include <GEL/CGLA/Vec.h>
#include <GEL/Geometry/QEM.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>
#include <GEL/HMesh/quadric_simplify.h>
#include <GEL/HMesh/smooth.h>


namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    using namespace Geometry;

    namespace
    {
        /** The simpliciation record contains information about a potential
         edge contraction */
        struct SimplifyRec
        {
            Vec3d opt_pos;  // optimal vertex position
            HalfEdgeID h;   // Index (into HalfEdgeRec vector) of edge
            // we want to contract
            float err;      // Error associated with contraction
            int time_stamp; // Time stamp (see comment on HalfEdgeRec)
            SimplifyRec(const Vec3d& _opt_pos, HalfEdgeID _h, float _err, int _time_stamp):
            opt_pos(_opt_pos), h(_h), err(_err), time_stamp(_time_stamp) {}
        };
        
        bool operator<(const SimplifyRec& s1, const SimplifyRec& s2)
        {
            return s1.err > s2.err;
        }
        
        /** This test is inspired by Garland's Ph.D. thesis. We try
         to detect whether flipped triangles will occur by sort of
         ensuring that the new vertex is in the hull of the one rings
         of the vertices at either end of the edge being contracted
         
         I also had an additional check intended to ensure that poor valencies
         would not be introduced, but it seemed to be unnecessary.
         
         */
        bool check_consistency(const Manifold& m, HalfEdgeID h, const Vec3d& opt_pos)
        {
            Walker w = m.walker(h);
            
            VertexID v0 = w.vertex();
            VertexID v1 = w.opp().vertex();
            Vec3d p0(m.pos(v0));
            
            for(Walker w = m.walker(v0); !w.full_circle(); w = w.circulate_vertex_cw()) {
                if(w.vertex()!= v1 && w.next().vertex() != v1){
                    Vec3d pa(m.pos(w.vertex()));
                    Vec3d pb(m.pos(w.next().vertex()));
                    
                    Vec3d dir = normalize(pb - pa);
                    
                    Vec3d n = p0 - pa;
                    n = n - dir * dot(dir,n);
                    
                    if(dot(n,opt_pos - pa) <= 0)
                        return false;
                }
            }
            return true;
        }
        
        
        /** This class contains the innards of the Garland-Heckbert simplification. The constructor
         makes a queue and the reduce function performa a number of reductions.
         */
        class SimplifyQueue {
            priority_queue<SimplifyRec> sim_queue;
            VertexAttributeVector<QEM> qem_vec;
            HalfEdgeAttributeVector<int> time_stamp;
            double singular_thresh;
            Manifold* m_ptr;
            
            SimplifyRec create_simplify_rec(HalfEdgeID h);
            
        public:
            SimplifyQueue(Manifold& m, double _singular_thresh);
            
            void reduce(long int max_work, double err_thresh);
        };
        
        SimplifyQueue::SimplifyQueue(Manifold& m, double _singular_thresh):
        m_ptr(&m),
        singular_thresh(_singular_thresh) {
            // For all vertices, compute quadric and store in qem_vec
            const auto processor_count = std::thread::hardware_concurrency();
            qem_vec = VertexAttributeVector<QEM>(m.allocated_vertices(), QEM());
            vector<thread> thread_vec(processor_count);
            for(int i=0;i<processor_count; ++i)
                thread_vec[i] = thread([&](int i){
                    for(auto v: m.vertices())
                        if (v.get_index()% processor_count == i) {
                            Vec3d p(m.pos(v));
                            Vec3d vn(normal(m, v));
                            QEM q;
                            for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
                                FaceID f = w.face();
                                if(f != InvalidFaceID){
                                    Vec3d n(normal(m, f));
                                    double a = area(m, f);
                                    q += QEM(p, n, a / 3.0);
                                }
                                if ((f == InvalidFaceID || w.opp().face() == InvalidFaceID ) && sqr_length(vn) > 0.0){
                                    Vec3d edge = Vec3d(m.pos(w.vertex())) - p;
                                    double edge_len = sqr_length(edge);
                                    if(edge_len > 0.0){
                                        Vec3d n = cross(vn, edge);
                                        q += QEM(p, n, 2*edge_len);
                                    }
                                }
                            }
                            qem_vec[v] = q;
                        }
                }, i);
            
            for(int i=0;i<processor_count; ++i)
                thread_vec[i].join();
            
            vector<SimplifyRec> recs;
            for(HalfEdgeID h: m.halfedges()) {
                time_stamp[h] = 0;
                if(h<m.walker(h).opp().halfedge())
                    recs.push_back(create_simplify_rec(h));
            }
            sim_queue = priority_queue<SimplifyRec>(begin(recs), end(recs));
        }
        
        SimplifyRec SimplifyQueue::create_simplify_rec(HalfEdgeID h)
        {
            Walker w = m_ptr->walker(h);
            time_stamp[h] +=1;
            
            VertexID hv = w.vertex();
            VertexID hov = w.opp().vertex();
            
            QEM q = qem_vec[hv];
            q += qem_vec[hov];
            
            Vec3d opt_origin = Vec3d(m_ptr->pos(hv) + m_ptr->pos(hov)) * 0.5;
            Vec3d opt_pos = q.opt_pos(singular_thresh, opt_origin);
            
            // Create SimplifyRec
            return SimplifyRec(opt_pos, h, q.error(opt_pos), time_stamp[h]);
        }
        
        
        void SimplifyQueue::reduce(long int max_work, double err_thresh)
        {
            int work = 0;
            while(!sim_queue.empty() && work < max_work){
                SimplifyRec simplify_rec = sim_queue.top();
                sim_queue.pop();
                
                if (simplify_rec.err > err_thresh)
                    return;
                
                HalfEdgeID h = simplify_rec.h;
                // First we check that the edge has not been removed and that it is
                // still the lower numbered halfedge in the pair.
                if (m_ptr->in_use(h) && h < m_ptr->walker(h).opp().halfedge()) {
                    
                    // Check the time stamp to verify that the simplification
                    // record is the newest.
                    if(time_stamp[h] == simplify_rec.time_stamp) {
                        Walker w = m_ptr->walker(h);
                        Walker wo = w.opp();
                        VertexID v = wo.vertex();
                        VertexID n = w.vertex();
                        
                        // Check the edge is, in fact, collapsible
                        if(precond_collapse_edge(*m_ptr, h)){
                            // If our consistency checks pass, we are relatively
                            // sure that the contraction does not lead to a face flip.
                            if(check_consistency(*m_ptr, h, simplify_rec.opt_pos) &&
                               check_consistency(*m_ptr, wo.halfedge(), simplify_rec.opt_pos)){
                                qem_vec[n] += qem_vec[v];
                                m_ptr->collapse_edge(h);
                                m_ptr->pos(n) = simplify_rec.opt_pos;
                                for(Walker w = m_ptr->walker(n); !w.full_circle(); w = w.circulate_vertex_cw())
                                    sim_queue.push(create_simplify_rec(w.hmin()));
                                work += 1;
                            }
                        }
                    }
                }
            }
        }
        
    } // end of anonymous namespace


    void quadric_simplify(Manifold& m, double keep_fraction, double singular_thresh, double _err_thresh)
    {
        int n = m.no_vertices();
        int max_work = max(0, int(n - keep_fraction * n));
        SimplifyQueue sq(m, singular_thresh);
        Vec3d c;
        float r;
        bsphere(m, c, r);
        double err_thresh = sqr(_err_thresh*r);
        sq.reduce(max_work, err_thresh);
    }
    
}
