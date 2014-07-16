/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../CGLA/Vec4f.h"
#include "Triangle.h"


using namespace std;
using namespace CGLA;

namespace Geometry
{

const float EPSILON = 1e-10f;

Triangle::Triangle(const CGLA::Vec3f& _v0, 
									 const CGLA::Vec3f& _v1, 
									 const CGLA::Vec3f& _v2,
									 
									 const CGLA::Vec3f& _vn0,
									 const CGLA::Vec3f& _vn1,
									 const CGLA::Vec3f& _vn2,
									 
									 const CGLA::Vec3f& _en0,
									 const CGLA::Vec3f& _en1,
									 const CGLA::Vec3f& _en2)
{
	vert[0]  =_v0;
	vert[1]  =_v1;
	vert[2]  =_v2;

	vert_norm[0] = _vn0;
	vert_norm[1] = _vn1;
	vert_norm[2] = _vn2;

	edge_norm[0] = _en0;
	edge_norm[1] = _en1;
	edge_norm[2] = _en2;

	face_norm = normalize(cross(vert[1]-vert[0], vert[2]-vert[0]));
	for(int i=0;i<3;++i)
		{
			int j= (i+1)%3;
			edge[i] = vert[j]-vert[i];
			tri_plane_edge_norm[i] = cross(face_norm, edge[i]);
			edge_len[i] = edge[i].length();
		}
}


// Moellers method
bool Triangle::intersect(const CGLA::Vec3f& orig,
												 const CGLA::Vec3f& dir, float&t) const
{
	Vec3f tvec, pvec, qvec;
	float det,inv_det;

   /* begin calculating determinant - also used to calculate U parameter */
   pvec = cross(dir, -edge[2]);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = dot(edge[0], pvec);

   if (det > -EPSILON && det < EPSILON)
     return 0;
   inv_det = 1.0 / det;

   /* calculate distance from v0 to ray origin */
   tvec =  orig - vert[0];

   /* calculate U parameter and test bounds */
   float u = dot(tvec, pvec) * inv_det;
   if (u < 0.0 || u > 1.0)
     return false;

   /* prepare to test V parameter */
   qvec = cross(tvec, edge[0]);

   /* calculate V parameter and test bounds */
   float v = dot(dir, qvec) * inv_det;
   if (v < 0.0 || u + v > 1.0)
     return false;

   /* calculate t, ray intersects triangle */
   t = dot(-edge[2], qvec) * inv_det;

   return true;
}




bool Triangle::signed_distance(const Vec3f& p, 
															 float& sq_dist, float& sgn) const
{
	int vertex_scores[3] = {0,0,0};
	Vec3f closest_pnt, normal;
	int idx_0;

	// Loop over all three edges.
	for(idx_0=0; idx_0<3; ++idx_0)
		{
			const int idx_1 = (idx_0+1) % 3;
			const Vec3f dir_3d = edge[idx_0]/edge_len[idx_0];
			const float t = dot(p - vert[idx_0], dir_3d);
			if(t <= 0)
				{
					++vertex_scores[idx_0];
					if(vertex_scores[idx_0] == 2)
						{
							closest_pnt = vert[idx_0];
							normal = vert_norm[idx_0];
							break;
						}
				}
			else if(t >= edge_len[idx_0])
				{
					++vertex_scores[idx_1];
					if(vertex_scores[idx_1] == 2)
						{
							closest_pnt = vert[idx_1];
							normal = vert_norm[idx_1];
							break;
						}
				}
			else if(dot(tri_plane_edge_norm[idx_0], p-vert[idx_0]) <=0)
				{
					closest_pnt=vert[idx_0]+t*dir_3d;
					normal = edge_norm[idx_0];
					break;
				}
		}
	if(idx_0 == 3)
		{
			closest_pnt = p - face_norm*(dot(p-vert[0],face_norm));
			normal = face_norm;
		}
	sq_dist = sqr_length(p-closest_pnt);
	
	// Compute dot product with angle weighted normal, and
	// assign the sign based on the result.
	if(dot(normal, p-closest_pnt) >=0)
		sgn = 1.0f;
	else
		sgn = -1.0f;

	return true;
}

}
