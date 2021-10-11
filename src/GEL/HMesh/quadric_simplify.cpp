/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <queue>
#include <iostream>
#include <thread>
#include <algorithm>
#include <random>

#include <GEL/CGLA/Vec3d.h>
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

        /* The simpliciation record contains information about a potential
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
            // record).
        };

        bool operator<(const SimplifyRec& s1, const SimplifyRec& s2)
        {
            return s1.err > s2.err;
        }

        class QuadricSimplifier
        {
            typedef priority_queue<SimplifyRec> SimplifyQueue;
            typedef VertexAttributeVector<QEM> QEMVec;
            typedef HalfEdgeAttributeVector<int> HalfEdgeTimeStamp;
            
            Manifold& m;
            HalfEdgeTimeStamp halfedge_time_stamp;
            QEMVec qem_vec;
            SimplifyQueue sim_queue;
            double singular_thresh;
            bool choose_optimal_positions;
            
            /* Compute the error associated with contraction of he and the
             optimal position of resulting vertex. */
            SimplifyRec create_simplify_rec(HalfEdgeID h);
            
            /* Check whether the contraction is valid. See below for details*/
            bool check_consistency(HalfEdgeID h, const Vec3d& opt_pos);
            
            /* Update time stamps for all halfedges in one ring of vi */
            void update_onering_timestamp(VertexID v);
            
            /* Perform a collapse - if conditions are met. Returns 1 or 0
             accordingly. */
            int collapse(SimplifyRec& simplify_rec);
            
        public:
            
            /* Create a simplifier for a manifold */
            QuadricSimplifier(Manifold& _m, double _singular_thresh, bool _choose_optimal_positions):
            m(_m),
            singular_thresh(_singular_thresh),
            choose_optimal_positions(_choose_optimal_positions)
            {}
            
            /* Simplify doing at most max_work contractions */
            void reduce(long int max_work);
        };

        bool QuadricSimplifier::check_consistency(HalfEdgeID h, const Vec3d& opt_pos)
        {
            Walker w = m.walker(h);
            
            
            VertexID v0 = w.vertex();
            VertexID v1 = w.opp().vertex();
            Vec3d p0(m.pos(v0));
            
            /* This test is inspired by Garland's Ph.D. thesis. We try
             to detect whether flipped triangles will occur by sort of
             ensuring that the new vertex is in the hull of the one rings
             of the vertices at either end of the edge being contracted
             
             I also had an additional check intended to ensure that poor valencies
             would not be introduced, but it seemed to be unnecessary.
             
             */
            
            for(Walker w = m.walker(v0); !w.full_circle(); w = w.circulate_vertex_cw()){
                //ConstHalfEdgeHandle h = vc.halfedge();
                
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


        SimplifyRec QuadricSimplifier::create_simplify_rec(HalfEdgeID h)
        {
            Walker w = m.walker(h);
            halfedge_time_stamp[h] +=1;
            
            VertexID hv = w.vertex();
            VertexID hov = w.opp().vertex();
            
            QEM q = qem_vec[hv];
            q += qem_vec[hov];
            
            Vec3d opt_pos(0);
            Vec3d opt_origin = Vec3d(m.pos(hv) + m.pos(hov)) * 0.5;
            if(choose_optimal_positions)
                opt_pos = q.opt_pos(singular_thresh,opt_origin);
            else
                opt_pos = Vec3d(m.pos(hv));
            
            float err = q.error(opt_pos);
            
            // Create SimplifyRec
            return SimplifyRec(opt_pos, h, err, halfedge_time_stamp[h]);
        }


        void QuadricSimplifier::update_onering_timestamp(VertexID v)
        {
            // For all emanating edges h
            for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()) {
                HalfEdgeID h = w.hmin();
                auto rec = create_simplify_rec(h);
                sim_queue.push(rec);
            }
        }

        int QuadricSimplifier::collapse(SimplifyRec& simplify_rec)
        {
            HalfEdgeID h = simplify_rec.h;
            if (m.in_use(h) &&
                h < m.walker(h).opp().halfedge()) {
                
                // Check the time stamp to verify that the simplification
                // record is the newest. If the halfedge has been removed
                // the time stamp is -1 and the comparison will also fail.
                if(halfedge_time_stamp[h] == simplify_rec.time_stamp) {
                    Walker w = m.walker(h);
                    Walker wo = w.opp();
                    VertexID v = wo.vertex();
                    VertexID n = w.vertex();
                    
                    // If the edge is, in fact, collapsible
                    if(precond_collapse_edge(m, h)){
                        // If our consistency checks pass, we are relatively
                        // sure that the contraction does not lead to a face flip.
                        if(check_consistency(h, simplify_rec.opt_pos) && check_consistency(wo.halfedge(), simplify_rec.opt_pos)){
                            qem_vec[n] += qem_vec[v];
                            m.collapse_edge(h);
                            m.pos(n) = simplify_rec.opt_pos;
                            update_onering_timestamp(n);
                            return 1;
                        }
                    }
                }
            }
            
            return 0;
        }

        void QuadricSimplifier::reduce(long int max_work)
        {
            cout << "Computing QEMs" << endl;

            // For all vertices, compute quadric and store in qem_vec
            const auto processor_count = std::thread::hardware_concurrency();
            qem_vec = QEMVec(m.allocated_vertices(), QEM());
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
            
            cout << "Pushing initial halfedges" << endl;
            
            vector<SimplifyRec> recs;
            for(HalfEdgeID h: m.halfedges()) {
                halfedge_time_stamp[h] = 0;
                if(h<m.walker(h).opp().halfedge())
                    recs.push_back(create_simplify_rec(h));
            }
            sim_queue = SimplifyQueue(begin(recs), end(recs));
            
            cout << "emptying queue" << endl;
            
            int work = 0;
            while(!sim_queue.empty() && work < max_work){
                SimplifyRec simplify_record = sim_queue.top();
                sim_queue.pop();
                work += collapse(simplify_record);
            }
        }
    }


    void quadric_simplify(Manifold& m, double keep_fraction, double singular_thresh, bool choose_optimal_positions)
    {
        int n = m.no_vertices();
        int max_work = max(0, int(n - keep_fraction * n));
        QuadricSimplifier qsim(m, singular_thresh, choose_optimal_positions);
        qsim.reduce(max_work);
    }
    
}
