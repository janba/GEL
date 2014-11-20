/*
 *  harmonics.cpp
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 01/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#include "harmonics.h"

#if USE_SPARSE_MATRIX
#include <arlgsym.h>
#endif

#include <CGLA/Vec3d.h>
#include <CGLA/Mat2x2d.h>
#include <LinAlg/Matrix.h>
#include <LinAlg/Vector.h>
#include <LinAlg/LapackFunc.h>

#include <GL/glew.h>
#include <GLGraphics/glsl_shader.h>

#include <HMesh/mesh_optimization.h>
#include <HMesh/curvature.h>
#include <HMesh/triangulate.h>
#include <HMesh/load.h>
#include <HMesh/x3d_save.h>

using namespace CGLA;
using namespace std;
using namespace HMesh;
using namespace Geometry;
using namespace LinAlg;


void Harmonics::make_laplace_operator()
{
	Q.Resize(mani.no_vertices(), mani.no_vertices());
	
	for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		if(!boundary(mani, *v)){
			int i = vtouched[*v];
			double area_i = mixed_area(mani, *v);
			Vec3d vertex(mani.pos(*v));
			Vec3d curv_normal(0);
			double a_sum = 0;
            for(Walker wv = mani.walker(*v); !wv.full_circle(); wv = wv.circulate_vertex_cw())
			{
				int j = vtouched[wv.vertex()];
				double area_j = mixed_area(mani, wv.vertex());
				
				Vec3d nbr(mani.pos(wv.vertex()));
				Vec3d left(mani.pos(wv.next().vertex()));
				Vec3d right(mani.pos(wv.opp().prev().opp().vertex()));
				
				double d_left = dot(normalize(nbr-left),
									normalize(vertex-left));
				double d_right = dot(normalize(nbr-right),
									 normalize(vertex-right));
				double a_left  = acos(min(1.0, max(-1.0, d_left)));
				double a_right = acos(min(1.0, max(-1.0, d_right)));
				
				double w = 1.0/tan(a_left) + 1.0/tan(a_right);
				
				Q[i][j] = -w/sqrt(area_i*area_j);
				//Q[i][j] = -1;
				a_sum += Q[i][j];
			}
			Q[i][i] = -a_sum;
		}
	EigenSolutionsSym(Q,V);
}

vector<Vec3d> Harmonics::analyze_signal(const VertexAttributeVector<Vec3d>& sig) {
    vector<Vec3d> sig_proj(maximum_eigenvalue, Vec3d(0));
    
    for(int es=0; es<maximum_eigenvalue; ++es)
		for(VertexID v: mani.vertices())
			sig_proj[es] +=  sig[v] * Q[es][vtouched[v]];
    
    return sig_proj;
}

VertexAttributeVector<Vec3d> Harmonics::reconstruct_signal(const vector<Vec3d>& sig_proj, int max_es) {
    VertexAttributeVector<Vec3d> sig(mani.allocated_vertices(), Vec3d(0));
    for(int es=0; es<max_es; ++es)
        for(VertexID v: mani.vertices())
            sig[v] += sig_proj[es] * Q[es][vtouched[v]];
    return sig;
}



Harmonics::Harmonics(HMesh::Manifold& _mani):mani(_mani), vtouched(_mani.allocated_vertices(), 0)
{
  
    int i = 0;
    for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v, ++i)
        vtouched[*v] = i;
	maximum_eigenvalue = mani.no_vertices()-1;
	make_laplace_operator();
	
	proj.resize(maximum_eigenvalue);
	max_eig_values.resize(maximum_eigenvalue, 1e-10f);
	
	cout << "Projection magnitude" << endl;
	for(int es=0; es<maximum_eigenvalue; ++es)
	{
		proj[es] = Vec3d(0.0);
		for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		{
			proj[es] +=  Vec3d(mani.pos(*v)) * Q[es][vtouched[*v]];
			max_eig_values[es] = max(max_eig_values[es], static_cast<float>(abs(Q[es][vtouched[*v]])));
		}
	}
}

void Harmonics::add_frequency(int es, float scale)
{
    if(es<maximum_eigenvalue) {
//        double sauce =  pow(max(0.0,1.0-.10*V[es]/V[V.Length()-1]),scale);
//        cout <<sauce << endl;
		for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v){
			Vec3d p = Vec3d(proj[es]);
			double Qval = Q[es][vtouched[*v]];
			mani.pos(*v) += p * Qval * scale;
		}
    }
}

void Harmonics::reset_shape()
{
	for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		mani.pos(*v) = Vec3d(0);
}
void Harmonics::partial_reconstruct(int E0, int E1, float scale)
{
	for(int es=E0;es<=E1;++es)
		add_frequency(es, scale);
}


double Harmonics::compute_adf(HMesh::VertexAttributeVector<double>& F, double t, double fiedler_boost)
{
	double F_max = 0;
    for(VertexID vid : mani.vertices())
        F[vid]=0;
    
	for(VertexID vid : mani.vertices()){
		for(int e = 1; e < V.Length(); ++e)
			F[vid] += sqr(Q[e][vtouched[vid]]) * exp(-t*V[e]/V[1]);
            if(fiedler_boost>0)
                for(int e = 1; e < V.Length(); ++e)
                    F[vid] += sqr(Q[e][vtouched[vid]]) * exp(-t) * fiedler_boost;
        
		F_max = max(F[vid], F_max);
	}
    return F_max;
}

double Harmonics::compute_esum(HMesh::VertexAttributeVector<double>& F, int e0, int e1)
{
	double F_max = 0;
    for(VertexID vid : mani.vertices())
        F[vid]=0;
    
	for(VertexID vid : mani.vertices()){
		for(int e = e0; e < e1; ++e)
			F[vid] += sqr(Q[e][vtouched[vid]]);
        
		F_max = max(F[vid], F_max);
	}
    return F_max;
}
