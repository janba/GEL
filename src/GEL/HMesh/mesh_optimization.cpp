/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */


#include "mesh_optimization.h"

#include <cfloat>
#include <queue>
#include <vector>

#include "../CGLA/Vec3d.h"
#include "../Geometry/Implicit.h"
//#include "Manifold.h"
#include "AttributeVector.h"
#include "triangulate.h"
#include "smooth.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    using namespace Geometry;
	
    // Small utility functions
    namespace 
    {
        class LineSeg
        {
            const Vec3d p0;
            const Vec3d p1;
            const float l;
            const CGLA::Vec3d dir;
        public:
			
            LineSeg(const Vec3d& _p0, const Vec3d& _p1): 
			p0(_p0), p1(_p1), l(length(p1-p0)), dir((p1-p0)/l) {}
			
			bool inseg(const Vec3d& p) const
			{
				double t = dot(dir, p-p0);
				if(t<0)
					return false;
				if(t>l)
					return false;
				return true;
			}
        };
		
        Vec3d compute_normal(Vec3d* v)
        {
            Vec3d norm;
            for(int i = 0; i < 4; ++i)
            {
                norm[0] += (v[i][1]-v[(i+1)%4][1])*(v[i][2]+v[(i+1)%4][2]);
                norm[1] += (v[i][2]-v[(i+1)%4][2])*(v[i][0]+v[(i+1)%4][0]);
                norm[2] += (v[i][0]-v[(i+1)%4][0])*(v[i][1]+v[(i+1)%4][1]);
            }
            float l = norm.length();
            if(l>0.0f)
                norm /= l;
            return norm;
        }
		
        bool would_flip(const Manifold& m, HalfEdgeID h)
        {
            Walker w = m.walker(h);
			
            VertexID hv = w.vertex();
            VertexID hov = w.opp().vertex();
            VertexID hnv = w.next().vertex();
            VertexID honv = w.opp().next().vertex();
			
            Vec3d v[4];
            v[0] = Vec3d(m.pos(hv));
            v[1] = Vec3d(m.pos(hov));
            v[2] = Vec3d(m.pos(hnv));
            v[3] = Vec3d(m.pos(honv));
			
            Vec3d dir = compute_normal(v);
			
            Vec3d n1a = cross(v[3]-v[0], v[2]-v[0]);
            Vec3d n2a = cross(v[2]-v[1], v[3]-v[1]);
			
            if(dot(normalize(n1a), dir) < 0)
                return true;
            if(dot(normalize(n2a), dir) < 0)
                return true;
            return false;
        }
    }
	
	
    double ValencyEnergy::delta_energy(const Manifold& m, HalfEdgeID h) const
    {
        Walker w = m.walker(h);
        
        VertexID v1 = w.opp().vertex();
        VertexID v2 = w.vertex();
        VertexID vo1 = w.next().vertex();
        VertexID vo2 = w.opp().next().vertex();
        
        int val1  = valency(m, v1);
        int val2  = valency(m, v2);
        int valo1 = valency(m, vo1);
        int valo2 = valency(m, vo2);
        
        // The optimal valency is four for a boundary vertex
        // and six elsewhere.
        int t1 = boundary(m, v1) ? 4 : 6;
        int t2 = boundary(m, v2) ? 4 : 6;
        int to1 = boundary(m, vo1) ? 4 : 6;
        int to2 = boundary(m, vo2) ? 4 : 6;
        
        int before = 
        sqr(val1-t1)+sqr(val2-t2)+
        sqr(valo1-to1)+sqr(valo2-to2);
        int after = 
        sqr(valo1+1-to1)+sqr(val1-1-t1)+
        sqr(val2-1-t2)+sqr(valo2+1-to2);
        
        return static_cast<double>(after-before);
    }
	
	
	class RandomEnergy: public EnergyFun
	{
	public:
		double delta_energy(const Manifold& m, HalfEdgeID he) const
		{
			return static_cast<double>(gel_rand()/static_cast<float>(GEL_RAND_MAX));
		}
	};
	
	
	double MinAngleEnergy::min_angle(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2) const
	{
		Vec3d a = normalize(v1-v0);
		Vec3d b = normalize(v2-v1);
		Vec3d c = normalize(v0-v2);
		
		return min(dot(a,-c), min(dot(b,-a), dot(c,-b)));
	}
	double MinAngleEnergy::delta_energy(const Manifold& m, HalfEdgeID h) const
	{
		Walker w = m.walker(h);
		
		VertexID hv = w.vertex();
		VertexID hnv = w.next().vertex();
		VertexID hov= w.opp().vertex();
		VertexID honv = w.opp().next().vertex();
		
		Vec3d v0(m.pos(hv));
		Vec3d v1(m.pos(hnv));
		Vec3d v2(m.pos(hov));
		Vec3d v3(m.pos(honv));
		
		Vec3d n1a = normalize(cross(v1-v0,v3-v0));
		Vec3d n2a = normalize(cross(v3-v2,v1-v2));
		
		if(dot(n1a, n2a) > thresh){
			double before = min(min_angle(v0,v1,v2), min_angle(v0,v2,v3));
			double after = min(min_angle(v0,v1,v3), min_angle(v1,v2,v3));
			return -(after-before);
		}
		return DBL_MAX;
	}
	
	
	void DihedralEnergy::compute_angles(const Manifold & m, HalfEdgeID h) const
	{
		Walker w = m.walker(h);
		
		VertexID hv = w.vertex();
		VertexID hov = w.opp().vertex();
		VertexID hnv = w.next().vertex();
		VertexID honv = w.opp().next().vertex();
		
		Vec3d va(m.pos(hv));
		Vec3d vb(m.pos(hov));
		Vec3d vc(m.pos(hnv));
		Vec3d vd(m.pos(honv));
		
		FaceID fa = w.next().opp().face();
		FaceID fb = w.next().next().opp().face();
		FaceID fc = w.opp().next().opp().face();
		FaceID fd = w.opp().next().next().opp().face();
		
		Vec3d n1 = normalize(cross(vc-va, vb-va));
		Vec3d n2 = normalize(cross(vb-va, vd-va));
		
		Vec3d na = fa == InvalidFaceID ? Vec3d(0) : Vec3d(normal(m, fa));
		Vec3d nb = fb == InvalidFaceID ? Vec3d(0) : Vec3d(normal(m, fb));
		Vec3d nc = fc == InvalidFaceID ? Vec3d(0) : Vec3d(normal(m, fc));
		Vec3d nd = fd == InvalidFaceID ? Vec3d(0) : Vec3d(normal(m, fd));
		
		
		Vec3d fn1 = normalize(cross(vb-vc, vd-vc));
		Vec3d fn2 = normalize(cross(vd-vc, va-vc));
		
		ab_12 = cos_ang(n1,n2);
		ab_a1 = cos_ang(na,n1);
		ab_b1 = cos_ang(nb,n1);
		ab_2c = cos_ang(n2,nc);
		ab_2d = cos_ang(n2,nd);
		
		aa_12 = cos_ang(fn1,fn2);
		aa_b1 = cos_ang(nb,fn1);
		aa_c1 = cos_ang(nc, fn1);
		aa_2a = cos_ang(fn2, na);
		aa_2d = cos_ang(fn2,nd);
	}
	
	double DihedralEnergy::energy(const Manifold& m, HalfEdgeID h) const
	{
		Walker w = m.walker(h);
		
		FaceID hf = w.face();
		FaceID hof = w.opp().face();
		
		double a = cos_ang(Vec3d(normal(m, hf)), Vec3d(normal(m, hof)));
		
		VertexID hv = w.vertex();
		VertexID hov = w.opp().vertex();
		
		Vec3d va(m.pos(hv));
		Vec3d vb(m.pos(hov));
		
		if(use_alpha)
			return edge_alpha_energy(va,vb,a);
		
		return edge_c_energy(va,vb,a);
	}
	
	
	double DihedralEnergy::delta_energy(const Manifold& m, HalfEdgeID h) const
	{
		compute_angles(m, h);
		
		Walker w = m.walker(h);
		
		VertexID hv = w.vertex();
		VertexID hov = w.opp().vertex();
		VertexID hnv = w.next().vertex();
		VertexID honv = w.opp().next().vertex();
		
		Vec3d va(m.pos(hv));
		Vec3d vb(m.pos(hov));
		Vec3d vc(m.pos(hnv));
		Vec3d vd(m.pos(honv));
		
		if(use_alpha){
			double before = 
			edge_alpha_energy(va,vb,ab_12)
			+edge_alpha_energy(va,vc,ab_a1)
			+edge_alpha_energy(vc,vb,ab_b1)
			+edge_alpha_energy(vd,vb,ab_2c)
			+edge_alpha_energy(vd,va,ab_2d);
			
			double after = 
			edge_alpha_energy(vd,vc,aa_12)
			+edge_alpha_energy(vb,vc,aa_b1)
			+edge_alpha_energy(vd,vb,aa_c1)
			+edge_alpha_energy(va,vc,aa_2a)
			+edge_alpha_energy(vd,va,aa_2d);
			
			return (after-before);
		}
		double before = 
		edge_c_energy(va,vb,ab_12)
		+edge_c_energy(va,vc,ab_a1)
		+edge_c_energy(vc,vb,ab_b1)
		+edge_c_energy(vd,vb,ab_2c)
		+edge_c_energy(vd,va,ab_2d);
		
		double after = 
		edge_c_energy(vd,vc,aa_12)
		+edge_c_energy(vb,vc,aa_b1)
		+edge_c_energy(vd,vb,aa_c1)
		+edge_c_energy(va,vc,aa_2a)
		+edge_c_energy(vd,va,aa_2d);
		
		return after-before;
	}
	
	
	
	double CurvatureEnergy::abs_mean_curv(const Vec3d& v, const vector<Vec3d>& ring) const
	{
        const size_t N = ring.size();

        double H = 0;

        Vec3d vnim1 = ring[N-1] - v;
        Vec3d vni   = ring[0] - v;
        Vec3d vnip1 = ring[1] - v;
        Vec3d Nm;
        Vec3d Np = normalize(cross(vni, vnim1));;
        for(size_t i = 0; i < N; ++i){
            Nm = Np;
            Np = normalize(cross(vnip1, vni));

            double beta = acos(max(-1.0, min(1.0, dot(Nm, Np))));
            H += vni.length() * beta;

            vni = vnip1;
            vnip1 = ring[(i+2)%N] - v;
        }

        return H/4;
	}
	
	double CurvatureEnergy::delta_energy(const Manifold& m, HalfEdgeID h) const
	{
        Walker w = m.walker(h);

        VertexID va = w.vertex();
        VertexID vb = w.opp().vertex();
        VertexID vc = w.next().vertex();
        VertexID vd = w.opp().next().vertex();

        Vec3d va_pos(m.pos(va));
        Vec3d vb_pos(m.pos(vb));
        Vec3d vc_pos(m.pos(vc));
        Vec3d vd_pos(m.pos(vd));

        for(Walker wv = m.walker(va); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
            VertexID v = wv.vertex();
            Vec3d pos(m.pos(v));

            va_ring_bef.push_back(pos);
            if(v != vb)
                va_ring_aft.push_back(pos);
        }
        for(Walker wv = m.walker(vb); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
            VertexID v = wv.vertex();
            Vec3d pos(m.pos(v));

            vb_ring_bef.push_back(pos);
            if(v != va)
                vb_ring_aft.push_back(pos);
        }
        for(Walker wv = m.walker(vc); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
            VertexID v = wv.vertex();
            Vec3d pos(m.pos(v));

            vc_ring_bef.push_back(pos);
            vc_ring_aft.push_back(pos);
            if(v == va)
                vc_ring_aft.push_back(vd_pos);
        }
        for(Walker wv = m.walker(vd); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
            VertexID v = wv.vertex();
            Vec3d pos(m.pos(v));

            vd_ring_bef.push_back(pos);
            vd_ring_aft.push_back(pos);
            if(v == vb)
                vd_ring_aft.push_back(vc_pos);
        }
        double before =
                abs_mean_curv(va_pos, va_ring_bef) +
                        abs_mean_curv(vb_pos, vb_ring_bef) +
                        abs_mean_curv(vc_pos, vc_ring_bef) +
                        abs_mean_curv(vd_pos, vd_ring_bef);

        double after =
                abs_mean_curv(va_pos, va_ring_aft) +
                        abs_mean_curv(vb_pos, vb_ring_aft) +
                        abs_mean_curv(vc_pos, vc_ring_aft) +
                        abs_mean_curv(vd_pos, vd_ring_aft);

        va_ring_bef.clear();
        va_ring_aft.clear();
        vb_ring_bef.clear();
        vb_ring_aft.clear();
        vc_ring_bef.clear();
        vc_ring_aft.clear();
        vd_ring_bef.clear();
        vd_ring_aft.clear();
		
		return after-before;
	}
	
	
	class GaussCurvatureEnergy: public EnergyFun
	{
		
		double gauss_curv(const Vec3d& v, const vector<Vec3d>& ring) const
		{
			const size_t N = ring.size();
			double asum=0.0f;
			double area_sum=0;
			for(size_t i = 0; i < N; ++i){
				const Vec3d& v1 = ring[i];
				const Vec3d& v2 = ring[(i+1)%N];
				Vec3d a = v1-v;
				Vec3d b = v2-v;
				asum += acos(max(-1.0, min(1.0, dot(a,b)/(length(a)*length(b)))));
				area_sum += 0.5 * length(cross(a,b));
			}
			return 3*abs(2 * M_PI - asum)/area_sum;
		}
		
	public:
		double delta_energy(const Manifold& m, HalfEdgeID h) const
		{
			Walker w = m.walker(h);
			
			VertexID va = w.vertex();
			VertexID vb = w.opp().vertex();
			VertexID vc = w.next().vertex();
			VertexID vd = w.opp().next().vertex();
			
			Vec3d va_pos(m.pos(va));
			Vec3d vb_pos(m.pos(vb));
			Vec3d vc_pos(m.pos(vc));
			Vec3d vd_pos(m.pos(vd));
			
			vector<Vec3d> va_ring_bef;
			vector<Vec3d> va_ring_aft;
			vector<Vec3d> vb_ring_bef;
			vector<Vec3d> vb_ring_aft;
			vector<Vec3d> vc_ring_bef;
			vector<Vec3d> vc_ring_aft;
			vector<Vec3d> vd_ring_bef;
			vector<Vec3d> vd_ring_aft;
			
			for(Walker wv = m.walker(va); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
				VertexID v = wv.vertex();
				Vec3d pos(m.pos(v));
				
				va_ring_bef.push_back(pos);
				if(v != vb)
					va_ring_aft.push_back(pos);
			}
			for(Walker wv = m.walker(vb); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
				VertexID v = wv.vertex();
				Vec3d pos(m.pos(v));
				
				vb_ring_bef.push_back(pos);
				if(v != va)
					vb_ring_aft.push_back(pos);
			}
			for(Walker wv = m.walker(vc); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
				VertexID v = wv.vertex();
				Vec3d pos(m.pos(v));
				
				vc_ring_bef.push_back(pos);
				vc_ring_aft.push_back(pos);
				if(v == va)
					vc_ring_aft.push_back(pos);
			}
			for(Walker wv = m.walker(vd); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
				VertexID v = wv.vertex();
				Vec3d pos(m.pos(v));
				
				vd_ring_bef.push_back(pos);
				vd_ring_aft.push_back(pos);
				if(v == vb)
					vd_ring_aft.push_back(pos);
			}
			double before =
			gauss_curv(va_pos,va_ring_bef) +
			gauss_curv(vb_pos,vb_ring_bef) +
			gauss_curv(vc_pos,vc_ring_bef) +
			gauss_curv(vd_pos,vd_ring_bef);
			
			double after =
			gauss_curv(va_pos,va_ring_aft) +
			gauss_curv(vb_pos,vb_ring_aft) +
			gauss_curv(vc_pos,vc_ring_aft) +
			gauss_curv(vd_pos,vd_ring_aft);
			
			return after-before;
		}
	};
	
	struct HalfEdgeCounter {
        int touched;
        bool isRemovedFromQueue;
    };

    struct PQElement
    {
        double pri;
        HalfEdgeID h;
        int time;

        //PQElement() {}
        PQElement(double _pri, HalfEdgeID _h, int _time):
                pri(_pri), h(_h), time(_time) {}
    };

    bool operator<(const PQElement & e0, const PQElement & e1)
    {
        return e0.pri > e1.pri;
    }

    void add_to_queue(const Manifold& m, HalfEdgeAttributeVector<HalfEdgeCounter>& counter, priority_queue<PQElement>& Q, HalfEdgeID h, const EnergyFun& efun, VertexAttributeVector<int> & flipCounter, int time)
    {
        if(boundary(m, h))
            return;

        Walker w = m.walker(h);

        // only consider one of the halfedges
        if (w.vertex() < w.opp().vertex()){
            h = w.opp().halfedge();
        }
        // if half edge already tested for queue in the current frame then skip
        if (counter[h].touched == time){
            return;
        }
        counter[h].isRemovedFromQueue = false;

        if(!precond_flip_edge(m, h))
            return;

        double energy = efun.delta_energy(m, h);
        counter[h].touched = time;

        const int avgValence = 6;
        if((energy < 0) && (flipCounter[w.vertex()] < avgValence)){
            Q.push(PQElement(energy, h, time));
        }
    }



	void add_one_ring_to_queue(const Manifold& m, HalfEdgeAttributeVector<HalfEdgeCounter>& touched, priority_queue<PQElement>& Q, VertexID v, const EnergyFun& efun, VertexAttributeVector<int> & flipCounter, int time)
	{
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
            add_to_queue(m, touched, Q, w.halfedge(), efun, flipCounter, time);
        }
	}
	
	
	void priority_queue_optimization(Manifold& m, const EnergyFun& efun)
	{
        HalfEdgeAttributeVector<HalfEdgeCounter> counter(m.allocated_halfedges(), HalfEdgeCounter{0, false});
        VertexAttributeVector<int> flipCounter(m.allocated_vertices(), 0);
        priority_queue<PQElement> Q;
		
		cout << "Building priority queue"<< endl;
        int time=1;
        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
            if(!counter[*h].touched) {
                add_to_queue(m, counter, Q, *h, efun, flipCounter, time);
            }
        }
		
		cout << "Emptying priority queue of size: " << Q.size() << " ";
		while(!Q.empty())
		{
			if(Q.size() % 1000 == 0)
				cout << ".";
			if(Q.size() % 10000 == 0)
				cout << Q.size();

            PQElement elem = Q.top();
            Q.pop();

            Walker w = m.walker(elem.h);

            if(counter[elem.h].isRemovedFromQueue) // if item already has been processed continue
                continue;

            counter[elem.h].isRemovedFromQueue = true;

            if(counter[elem.h].touched != elem.time) {
                if (efun.delta_energy(m, elem.h) >= 0) {
                    continue;
                }
            }
            if(!precond_flip_edge(m, elem.h))
                continue;

            flipCounter[w.vertex()]++;

            m.flip_edge(elem.h);

            add_one_ring_to_queue(m, counter, Q, w.vertex(), efun, flipCounter, time);
            add_one_ring_to_queue(m, counter, Q, w.next().vertex(), efun, flipCounter, time);
            add_one_ring_to_queue(m, counter, Q, w.opp().vertex(), efun, flipCounter, time);
            add_one_ring_to_queue(m, counter, Q, w.opp().next().vertex(), efun, flipCounter, time);

        }
		cout << endl;
	}
	
	void simulated_annealing_optimization(Manifold& m, const EnergyFun& efun, int max_iter)
	{
		gel_srand(0);
		int swaps;
		int iter = 0;
		double T = 1;
		
		double tmin=0;
		{
			for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
				if(boundary(m, *h)) 
					continue;  
				double e = efun.delta_energy(m, *h);
				tmin = std::min(e, tmin);
			}
		}
		if (tmin < 0.0) 
			T = -2*tmin;
		
		if(max_iter>0){ 
			do{
				cout << "Temperature : " << T << endl;
				vector<HalfEdgeID>  halfedges;
				for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
					if(boundary(m, *h))
						continue;
					halfedges.push_back(*h);
				}
				random_shuffle(halfedges.begin(), halfedges.end());
				swaps = 0;
				for(size_t i = 0; i < halfedges.size(); ++i){
					HalfEdgeID h = halfedges[i];
					DihedralEnergy dih_en;
					double dma = dih_en.min_angle(m, h);
					if(dma > -0.4){
						double delta = efun.delta_energy(m, h);
						if(delta < -1e-8){
							if(precond_flip_edge(m, h)){
								m.flip_edge(h);
								++swaps;
							}
						}
						else{
							delta = max(1e-8, delta);
							double prob = min(0.9999, exp(-delta/T));
							if(gel_rand()/double(GEL_RAND_MAX) < prob){
								if(precond_flip_edge(m, h)){
									m.flip_edge(h);
									++swaps;
								}
							}
						}
					}
				}
				cout << "Swaps = " << swaps << " T = " << T << endl;
				if(iter % 5 == 0 && iter > 0)
					T *= 0.9;
			}
			while(++iter < max_iter && swaps);
		}
		cout << "Iterations "  << iter << endl; 
		
	}
	
	
	void minimize_dihedral_angle(Manifold& m, 
								 int iter,
								 bool anneal, 
								 bool alpha,
								 double gamma)
	{
		DihedralEnergy energy_fun(gamma, alpha);
		if(anneal)
			simulated_annealing_optimization(m, energy_fun, iter);
		else
			priority_queue_optimization(m, energy_fun);
	}
	
	void randomize_mesh(Manifold& m, int max_iter)
	{
		RandomEnergy energy_fun;
		simulated_annealing_optimization(m, energy_fun, max_iter);
	}
	
	void minimize_curvature(Manifold& m, bool anneal)
	{
		CurvatureEnergy energy_fun;
		if(anneal)
			simulated_annealing_optimization(m, energy_fun);
		else
			priority_queue_optimization(m, energy_fun);
	}
	
	void minimize_gauss_curvature(Manifold& m, bool anneal)
	{
		GaussCurvatureEnergy energy_fun;
		if(anneal)
			simulated_annealing_optimization(m, energy_fun);
		else
			priority_queue_optimization(m, energy_fun);
	}
	
	void maximize_min_angle(Manifold& m, float thresh, bool anneal)
	{
		MinAngleEnergy energy_fun(thresh);
		if(anneal)
			simulated_annealing_optimization(m, energy_fun);
		else
			priority_queue_optimization(m, energy_fun);
	}
	
	void optimize_valency(Manifold& m, bool anneal)
	{
		ValencyEnergy energy_fun;
		if(anneal)
			simulated_annealing_optimization(m, energy_fun);
		else
			priority_queue_optimization(m, energy_fun);
	}
    
    
    void edge_equalize(Manifold& m, const Implicit& imp, int max_iter)
    {
        
        for(int iter=0;iter<max_iter;++iter)
        {
            float avg_edge_len=0;
            for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h)
                avg_edge_len += length(m, *h);
            avg_edge_len /= m.no_halfedges();
            
            TAL_smoothing(m,1,1);
            for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
                imp.push_to_surface(m.pos(*vid),0,avg_edge_len*0.5);
            
            vector<float> edge_lengths;
            int n=0;
            for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h,++n)
                edge_lengths.push_back(length(m, *h));
            sort(edge_lengths.begin(), edge_lengths.end());
            float med_len = edge_lengths[n/2];
            
            vector<HalfEdgeID> short_edges, long_edges;
            for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h)
                if(*h<m.walker(*h).opp().halfedge())
                {
                    float l = length(m,*h);
                    if(l> (4/3.0) * med_len)
                        long_edges.push_back(*h);
                    else if(l< (4/5.0) * med_len)
                        short_edges.push_back(*h);
                }
            for(int i=0;i<long_edges.size(); ++i)
                if(m.in_use(long_edges[i]))
                    m.split_edge(long_edges[i]);
            shortest_edge_triangulate(m);
            
            for(int i=0;i<short_edges.size(); ++i)
                if(m.in_use(short_edges[i]) && precond_collapse_edge(m, short_edges[i]))
                {
                    Walker w = m.walker(short_edges[i]);
                    Vec3d mid_point = 0.5*(m.pos(w.vertex())+m.pos(w.opp().vertex()));
                    bool illegal = false;
                    for(Walker wc=m.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
                        if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* med_len)
                            illegal = true;
                    for(Walker wc=m.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
                        if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* med_len)
                            illegal = true;
                    
                    if(!illegal)
                    {
                        if(boundary(m,w.vertex()) && !boundary(m,w.opp().vertex()))
                            m.collapse_edge(short_edges[i],false);
                        else if(!boundary(m,w.vertex()) && boundary(m,w.opp().vertex()));
                        else
                            m.collapse_edge(short_edges[i],true);
                    }
                }
            //        TAL_smoothing(m,.95,5);
            //        for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            //            imp.push_to_surface(m.pos(*vid));
            
            maximize_min_angle(m, 0.9);
        }
        m.cleanup();
        
    }
}
