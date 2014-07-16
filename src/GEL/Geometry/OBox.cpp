/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <cfloat>
#include "../CGLA/statistics.h"
#include "../CGLA/eigensolution.h"
#include "../CGLA/Mat4x4f.h"
#include "AABox.h"
#include "OBox.h"
#include "Triangle.h"

using namespace std;
using namespace CGLA;

namespace
{

	Mat3x3f compute_rotation(const vector<Geometry::Triangle>& invec)
	{
		const int N_tri = invec.size();

		Mat3x3f C;
		float a_H = 0;
		Vec3f m_H(0);
		
		for(int i=0;i<N_tri;++i)
			{
				const Geometry::Triangle& tri = invec[i];
				
				float a_k = tri.area();
				a_H += a_k;
				
				Vec3f m_k = tri.centre();
				m_H += a_k * m_k;
				
				Vec3f p_k = tri.get_v0();
				Vec3f q_k = tri.get_v1();
				Vec3f r_k = tri.get_v2();
				
				Mat3x3f M,P,Q,R;
				outer_product(m_k,m_k, M);
				outer_product(p_k,p_k, P);
				outer_product(q_k,q_k, Q);
				outer_product(r_k,r_k, R);
				
				C += (a_k/12.0f) * (9*M+P+Q+R);
			}
		m_H /= a_H;
		C /= a_H;
		Mat3x3f M;
		outer_product(m_H, m_H, M);
		C -= M;

		Mat3x3f Q,L;
		const int n_eig = power_eigensolution(C,Q,L,2);
		
		Vec3f X = normalize(Q[0]);
		Vec3f Y = normalize(Q[1]);
		Vec3f Z;
		
		float xy_ortho = fabs(dot(X,Y));
		if(n_eig == 2 && xy_ortho < 0.3)
			{
				if(xy_ortho > CGLA::TINY)
					Y = normalize(Y-X*dot(X,Y));
				Z = normalize(cross(X,Y));
			}
		else if(n_eig==1)
			{
				Y = Vec3f(X[2],X[0],X[1]);
				Y = normalize(Y-X*dot(X,Y));
				Z = normalize(cross(X,Y));
			}
		else
			{
				X=Vec3f(1,0,0);
				Y=Vec3f(0,1,0);
				Z=Vec3f(0,0,1);
			}
		return Mat3x3f(X,Y,Z);
		
	}



}

namespace Geometry
{

bool OBox::intersect(const CGLA::Vec3f& p, const CGLA::Vec3f& d) const 
{
	Vec3f pr = R * p;
	Vec3f dr = R * d;
	return aabox.intersect(pr, dr);
}


OBox OBox::box_triangle(const Triangle& t)
{
	Vec3f e0 = t.get_v1()-t.get_v0();
	Vec3f e1 = t.get_v2()-t.get_v1();
	Vec3f e2 = t.get_v0()-t.get_v2();

	Vec3f X,Y,Z;
	if(sqr_length(e0) > sqr_length(e1))
		{
			if(sqr_length(e0) > sqr_length(e2))
				{
					X = normalize(e0);
					Y = normalize(e1 - X * dot(X, e1));
				}
			else
				{
					X = normalize(e2);
					Y = normalize(e0 - X * dot(X, e0));
				}
		}
	else
		{
			if(sqr_length(e1) > sqr_length(e2))
				{
					X = normalize(e1);
					Y = normalize(e2 - X * dot(X, e2));
				}
			else
				{
					X = normalize(e2);
					Y = normalize(e0 - X * dot(X, e0));
				}
		}
	Z = cross(X,Y);
	
	const Mat3x3f Rot(X,Y,Z);

	Vec3f p0 = Rot * t.get_v0();
	Vec3f p1 = Rot * t.get_v1();
	Vec3f p2 = Rot * t.get_v2();
	Vec3f pmin = v_min(p0, v_min(p1, p2));
	Vec3f pmax = v_max(p0, v_max(p1, p2));
	
	Vec3f centre_close = v_max(pmin, v_min(pmax, Rot * t.get_centre()));
	return OBox(Rot, AABox(pmin, pmax, centre_close));
}


OBox OBox::box_and_split(const std::vector<Triangle>& invec,
													 std::vector<Triangle>& lvec,
													 std::vector<Triangle>& rvec)
{
	// Obtain the rotation matrix for the OBB
	const Mat3x3f Rot = compute_rotation(invec);
	const int N_tri = invec.size();
	const int N_pts = 3*N_tri;

	// Compute the rotated set of points and the extents of the point aligned 
	// BBox.
	vector<Vec3f> pts(N_pts);
	Vec3f tri_pmin(FLT_MAX), tri_pmax(-FLT_MAX);
	for(int i=0;i<N_tri;++i)
		{
			const Triangle& tri = invec[i];
			
			int offs = 3*i;
			pts[offs  ] = Rot*tri.get_v0();
			pts[offs+1] = Rot*tri.get_v1();
			pts[offs+2] = Rot*tri.get_v2();
			
			for(int j=0;j<3;++j)
				{
					tri_pmin = v_min(pts[offs+j], tri_pmin);
					tri_pmax = v_max(pts[offs+j], tri_pmax);
				}
		}

	// Find the point closest to the centre.
	const Vec3f centre = tri_pmin + 0.5f*(tri_pmax-tri_pmin);
	Vec3f centre_close;
	float min_dist = FLT_MAX;
	for(int i=0;i<N_pts;++i)
		{
			Vec3f v = pts[i];
			float sl = sqr_length(centre-v);
			if(sl < min_dist)
				{
					min_dist = sl;
					centre_close = v;
				}
		}

	// Partition the triangles
	const float thresh = centre[0];
	for(int i=0;i<N_tri;++i)
		{
			Vec3f p = Rot * invec[i].get_centre();
			if( p[0] > thresh)
				rvec.push_back(invec[i]);
			else
				lvec.push_back(invec[i]);
		}

	// If all triangles landed in one box, split them naively.
	if(lvec.empty() || rvec.empty())
		{
			lvec.clear();
			lvec.insert(lvec.end(),
									invec.begin(),
									invec.begin()+N_tri/2);
			rvec.clear();
			rvec.insert(rvec.end(),
									invec.begin()+N_tri/2,
									invec.end());
		}
		
	return OBox(Rot, AABox(tri_pmin, tri_pmax, centre_close));
}

}
