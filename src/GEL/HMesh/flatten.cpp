/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "flatten.h"

#include <string>
#include <fstream>
#include <vector>
#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
	
    void flatten(Manifold& m, WeightScheme ws)
    {
        HalfEdgeAttributeVector<double> edge_weights(m.allocated_halfedges(), 0);
		for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
		{
			for(Walker wv = m.walker(*f); !wv.full_circle(); wv = wv.circulate_face_ccw())
			{
				HalfEdgeID h = wv.halfedge();
				Vec3d p1(m.pos(wv.vertex()));
				Vec3d p2(m.pos(wv.next().vertex()));
				Vec3d p0(m.pos(wv.opp().vertex()));
				
				if(ws == FLOATER_W){
					double ang = acos(min(1.0, max(-1.0, dot(normalize(p1-p0), normalize(p2-p0)))));
					double ang_opp = acos(min(1.0, max(-1.0, dot(normalize(p2-p1), normalize(p0-p1)))));
					double l = (p1-p0).length();
					edge_weights[h]  += tan(ang/2) / l;
					edge_weights[wv.opp().halfedge()]  += tan(ang_opp/2) / l;
				}
				else if(ws == HARMONIC_W || ws == LSCM_W){
					double a = acos(min(1.0, max(-1.0, dot(normalize(p0-p2), normalize(p1-p2)))));
					double w = max(0.0000001,0.5/tan(a));
					edge_weights[h]  += w;
					edge_weights[wv.opp().halfedge()]  += w;
				}
				else{
					edge_weights[h]  = valency(m, wv.opp().vertex());
					edge_weights[wv.opp().halfedge()]  = valency(m, wv.vertex());
				}
			}
			
		}
		
		
        ofstream ofs("parametrized.obj");
		
        ofs << "mtllib parametrized.mtl\nusemtl mat\n" << endl;
		
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
            ofs << "v " << m.pos(*v)[0] << " " << m.pos(*v)[1] << " " << m.pos(*v)[2] << endl;
        ofs << endl;
		
		VertexAttributeVector<double> touched(m.allocated_vertices(), 0);
        VertexIDIterator v = m.vertices_begin();
        for(; v != m.vertices_end(); ++v){
            if(boundary(m, *v))
                break;
        }
        int n = 0;
        Walker bv = m.walker(*v);
        do{
            ++n;
            bv = bv.next();
        }
        while(bv.vertex() != *v);
		
        int i = 0;
        do{
			if(i==int(n*0.25) || i==int(n*0.75))
				touched[bv.vertex()]=1;
            double a = 2.0*M_PI*double(i)/n;
            m.pos(bv.vertex()) = Vec3d(cos(a), sin(a), 0);
            ++i;
            bv = bv.next();
        }
        while(bv.vertex() != *v);
		
        for(v = m.vertices_begin(); v != m.vertices_end(); ++v)
            if(!boundary(m, *v))
                m.pos(*v) = Vec3d(0.0);
		
        VertexAttributeVector<Vec3d> new_pos(m.no_vertices());
        for(int i = 0; i < 15000; ++i){
            for(v = m.vertices_begin(); v != m.vertices_end(); ++v){
				if(boundary(m, *v))
				{
					if(i>5000 && ws == LSCM_W && touched[*v] != 1)
					{
						Vec3d p_new(0);
						double w_sum = 1e-6;
						Vec3d grad_sum(0.0);
						Walker wv = m.walker(*v);
						for(;!wv.full_circle(); wv = wv.circulate_vertex_ccw())
                        {
							if(wv.face() != InvalidFaceID)
							{
								Vec3d p1(m.pos(wv.next().vertex()));
								Vec3d p0(m.pos(wv.vertex()));
								Vec3d area_grad = 0.5*(p1 - p0); 
								grad_sum[0] += -area_grad[1];
								grad_sum[1] += area_grad[0];
							}
                            double w = edge_weights[wv.halfedge()];
                            p_new += Vec3d(m.pos(wv.vertex()) * w);
                            w_sum += w;                            
                        }
						new_pos[*v] = ((p_new) - (grad_sum))/w_sum;	
					}
                    else
                        new_pos[*v] = m.pos(*v);
				}
				else
				{
					Vec3d p_new(0);
					double w_sum = 1e-6;
					for(Walker wv = m.walker(*v); !wv.full_circle(); wv = wv.circulate_vertex_ccw()) 
					{
						double w = edge_weights[wv.halfedge()];
						p_new += Vec3d(m.pos(wv.vertex()) * w);
						w_sum += w;
					}
                    new_pos[*v] = p_new/w_sum;
                }
            }
            for(v = m.vertices_begin(); v != m.vertices_end(); ++v)
                m.pos(*v) = new_pos[*v];
        }
		
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);
        i = 0;
        for(v = m.vertices_begin(); v != m.vertices_end(); ++v, ++i){
            ofs << "vt " << (0.5*m.pos(*v)[0]+0.5) << " " << (0.5*m.pos(*v)[1]+0.5)  << endl;
            vtouched[*v] = i;
        }
		
        ofs << endl;
		
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            ofs << "f ";
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_cw()){
                int idx = vtouched[w.vertex()] + 1;
                ofs << idx << "/" << idx <<" ";
            }
            ofs << endl;
        }
		
    }
}