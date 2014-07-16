/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "quadric_simplify.h"

#include <queue>
#include <iostream>
#include "../CGLA/Vec3d.h"
#include "../Geometry/QEM.h"

#include "Manifold.h"
#include "AttributeVector.h"
#include "smooth.h"


namespace HMesh
{
	using namespace std;
	using namespace CGLA;
    using namespace Geometry;
    
    namespace
    {
        /* We create a record for each halfedge where we can keep its time
         stamp. If the time stamp on the halfedge record is bigger than
         the stamp on the simplification record, we cannot use the 
         simplification record (see below). */
        struct HalfEdgeRec
        {
            HalfEdgeID h;
            int time_stamp;
            void halfedge_removed() {time_stamp = -1;}
            HalfEdgeRec(): time_stamp(0) {}
        };
        
        /* The simpliciation record contains information about a potential
         edge contraction */
        struct SimplifyRec
        {
            Vec3d opt_pos;  // optimal vertex position
            HalfEdgeID h;   // Index (into HalfEdgeRec vector) of edge
            // we want to contract
            float err;      // Error associated with contraction
            int time_stamp; // Time stamp (see comment on HalfEdgeRec)
            int visits;     // Visits (number of times we considered this 
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
            typedef HalfEdgeAttributeVector<HalfEdgeRec> HalfEdgeVec;
            typedef VertexAttributeVector<int> CollapseMask;
            
            Manifold& m;
            HalfEdgeVec halfedge_vec;
            QEMVec qem_vec;
            CollapseMask collapse_mask;
            SimplifyQueue sim_queue;
            double singular_thresh;
            bool choose_optimal_positions;
            
            /* Compute the error associated with contraction of he and the
             optimal position of resulting vertex. */
            void push_simplify_rec(HalfEdgeID h);
            
            /* Check whether the contraction is valid. See below for details*/
            bool check_consistency(HalfEdgeID h, const Vec3d& opt_pos);
            
            /* Update the time stamp of a halfedge. A halfedge and its opp edge
             may have different stamps. We choose a stamp that is greater
             than either and assign to both.*/
            void update_time_stamp(HalfEdgeID h)
            {
				Walker w = m.walker(h);
				HalfEdgeID ho = w.opp().halfedge();
                
                int time_stamp = std::max( halfedge_vec[h].time_stamp, halfedge_vec[ho].time_stamp);
                ++time_stamp;
                halfedge_vec[h].time_stamp = time_stamp;
                halfedge_vec[ho].time_stamp = time_stamp;
            }
            
            /* Update time stamps for all halfedges in one ring of vi */
            void update_onering_timestamp(VertexID v);
            
            /* Perform a collapse - if conditions are met. Returns 1 or 0 
             accordingly. */
            int collapse(SimplifyRec& simplify_rec);
            
        public:
            
            /* Create a simplifier for a manifold */
            QuadricSimplifier(Manifold& _m, VertexAttributeVector<int>& _collapse_mask,
                              double _singular_thresh, bool _choose_optimal_positions):
            m(_m), 
            halfedge_vec(_m.allocated_halfedges()), 
            qem_vec(_m.allocated_vertices()),
            collapse_mask(_collapse_mask),
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
        
        
        void QuadricSimplifier::push_simplify_rec(HalfEdgeID h)
        {
			Walker w = m.walker(h);
            if(collapse_mask[w.opp().vertex()] == 0)
            {
                
                update_time_stamp(h);
                
                VertexID hv = w.vertex();
                VertexID hov = w.opp().vertex();
                
                // Get QEM for both end points
                const QEM& Q1 = qem_vec[hv];
                const QEM& Q2 = qem_vec[hov];
                
                QEM q = Q1;
                q += Q2;
                
                float err;
                Vec3d opt_pos(0);
                Vec3d opt_origin = Vec3d(m.pos(hv) + m.pos(hov)) * 0.5;
                if(choose_optimal_positions)
                    opt_pos = q.opt_pos(singular_thresh,opt_origin);
                else
                    opt_pos = Vec3d(m.pos(hv));
                
                err = q.error(opt_pos);
                
                // Create SimplifyRec
                SimplifyRec simplify_rec;
                simplify_rec.opt_pos = opt_pos;
                simplify_rec.err = err;
                simplify_rec.h = h;
                simplify_rec.time_stamp = halfedge_vec[h].time_stamp;
                simplify_rec.visits = 0;
                // push it.
                sim_queue.push(simplify_rec);
            }
        }
        
        
        void QuadricSimplifier::update_onering_timestamp(VertexID v)
        {
            // For all emanating edges h
			for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw())
                push_simplify_rec(w.halfedge());
        }
        
        int QuadricSimplifier::collapse(SimplifyRec& simplify_rec)
        { 
            HalfEdgeID h = halfedge_vec[simplify_rec.h].h;
            
			Walker w = m.walker(h);
            
            // Check the time stamp to verify that the simplification 
            // record is the newest. If the halfedge has been removed
            // the time stamp is -1 and the comparison will also fail.
            if(halfedge_vec[h].time_stamp == simplify_rec.time_stamp){
				Walker wo = w.opp();
                VertexID v = wo.vertex();
                VertexID n = w.vertex();
                
                // If the edge is, in fact, collapsible
                if(precond_collapse_edge(m, h)){      
                    // If our consistency checks pass, we are relatively
                    // sure that the contraction does not lead to a face flip.
                    if(check_consistency(h, simplify_rec.opt_pos) && check_consistency(wo.halfedge(), simplify_rec.opt_pos)){
                        //cout << simplify_rec.err << " " << &(*he->vert) << endl;
                        // Get QEM for both end points
                        const QEM& Q1 = qem_vec[n];
                        const QEM& Q2 = qem_vec[v];
                        
                        // Compute Q_new = Q_1 + Q_2
                        QEM q = Q1;
                        q += Q2;
                        
                        // Mark all halfedges that will be removed as dead
                        halfedge_vec[w.halfedge()].halfedge_removed();
                        halfedge_vec[wo.halfedge()].halfedge_removed();
                        
                        if(w.next().next().next().halfedge() == h){
                            halfedge_vec[w.next().halfedge()].halfedge_removed();
                            halfedge_vec[w.next().next().halfedge()].halfedge_removed();
                        }
                        if(wo.next().next().next().halfedge() == wo.halfedge()){
                            halfedge_vec[wo.next().halfedge()].halfedge_removed();
                            halfedge_vec[wo.next().next().halfedge()].halfedge_removed();
                        }
                        
                        // Do collapse
                        m.collapse_edge(h);
                        m.pos(n) = simplify_rec.opt_pos;
                        qem_vec[n] = q;
                        
                        update_onering_timestamp(n);
                        return 1;
                    }
                }
                
                
                // If we are here, the collapse was not allowed. If we have
                // seen this simplify record less than 100 times, we try to
                // increase the error and store the record again. Maybe some
                // other contractions will make it more digestible later.
                if(simplify_rec.visits < 100){
                    simplify_rec.err *= 1.01f;
                    ++simplify_rec.visits;
                    sim_queue.push(simplify_rec);
                }
            }
            
            return 0;
        }
        
        void QuadricSimplifier::reduce(long int max_work)
        {
            // Set t = 0 for all halfedges
            for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
                halfedge_vec[*h].h = *h;
            }
            cout << "Computing quadrics" << endl;
            
            // For all vertices, compute quadric and store in qem_vec
            for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                Vec3d p(m.pos(*v));
                Vec3d vn(normal(m, *v));
                QEM q;
				for(Walker w = m.walker(*v); !w.full_circle(); w = w.circulate_vertex_cw()){
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
                qem_vec[*v] = q;
            }
            cout << "Pushing initial halfedges" << endl;
            
            for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
                if(halfedge_vec[*h].time_stamp == 0)
                    push_simplify_rec(*h);
            }
            cout << "Simplify";
            
            int work = 0;
            while(!sim_queue.empty() && work < max_work){
                SimplifyRec simplify_record = sim_queue.top();
                sim_queue.pop();
                
                work += 2*collapse(simplify_record);
                if((sim_queue.size() % 10000) == 0){
                    cout << ".";
                }
//                cout << "work = " << work << endl;
//                cout << "sim Q size = " << sim_queue.size() << endl;
            }
            cout << endl;
        }
    }
    
    void quadric_simplify(Manifold& m, double keep_fraction, double singular_thresh, bool choose_optimal_positions)
    {
        gel_srand(1210);
        long int F = m.no_faces();
        VertexAttributeVector<int> mask(m.no_faces(), 0);
        long int max_work = max(static_cast<long int>(0), F- static_cast<long int>(keep_fraction * F));
        QuadricSimplifier qsim(m, mask, singular_thresh, choose_optimal_positions);
        qsim.reduce(max_work);
    }
    
    void quadric_simplify(Manifold& m, VertexAttributeVector<int> mask, double keep_fraction, double singular_thresh, bool choose_optimal_positions)
    {
        gel_srand(1210);
        long int F = m.no_faces();
        long int max_work = keep_fraction == 0.0 ? INT_MAX : max(static_cast<long int>(0),
                                                                 F- static_cast<long int>(keep_fraction * F));
        QuadricSimplifier qsim(m, mask, singular_thresh, choose_optimal_positions);
        qsim.reduce(max_work);
    }
    
    
}
