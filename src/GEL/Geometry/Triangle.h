/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Triangle.h
 * @brief Triangle class for use with bbox hierarchies.
 */

#ifndef __GEOMETRY_TRIANGLE_H
#define __GEOMETRY_TRIANGLE_H

#include "../CGLA/Vec2f.h"
#include "../CGLA/Vec3f.h"
#include "../CGLA/Mat3x3f.h"

namespace Geometry
{

class Triangle
{
	CGLA::Vec3f vert[3];
	CGLA::Vec3f edge[3];

	CGLA::Vec3f vert_norm[3];
	CGLA::Vec3f edge_norm[3];
	CGLA::Vec3f face_norm;

	CGLA::Vec3f tri_plane_edge_norm[3];
	float edge_len[3];
	
 public:

	Triangle() {}
	
	Triangle(const CGLA::Vec3f& _v0, 
					 const CGLA::Vec3f& _v1, 
					 const CGLA::Vec3f& _v2,

					 const CGLA::Vec3f& _vn0,
					 const CGLA::Vec3f& _vn1,
					 const CGLA::Vec3f& _vn2,

					 const CGLA::Vec3f& _en0,
					 const CGLA::Vec3f& _en1,
					 const CGLA::Vec3f& _en2);

	bool intersect(const CGLA::Vec3f&, const CGLA::Vec3f&, float&) const;

	const CGLA::Vec3f get_pmin() const 
	{
		return v_min(vert[0],v_min(vert[1],vert[2]));
	}

	const CGLA::Vec3f get_pmax() const 
	{
		return v_max(vert[0],v_max(vert[1],vert[2]));
	}

	const CGLA::Vec3f centre() const
		{
			return (vert[0]+vert[1]+vert[2])/3.0f;
		}

	const float area() const 
		{
			return 0.5 * (cross(edge[0],-edge[2])).length();
		}

	const CGLA::Vec3f get_centre() const 
	{
		CGLA::Vec3f pmin = get_pmin();
		CGLA::Vec3f pmax = get_pmax();
		return (pmax-pmin)/2.0f+pmin;
	}

	bool signed_distance(const CGLA::Vec3f& p, float& sq_dist, float& sgn) const;

	CGLA::Vec3f get_v0() const {return vert[0];}
	CGLA::Vec3f get_v1() const {return vert[1];}
	CGLA::Vec3f get_v2() const {return vert[2];}

	const CGLA::Vec3f& get_edge(int i) const {return edge[i];}

	const CGLA::Vec3f& get_face_norm() const {return face_norm;}
	
};	

}
#endif
