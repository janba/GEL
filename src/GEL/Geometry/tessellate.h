/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file tessellate.h
 * @brief Tessellation of a parametric surface. Based on algorithm by Velho.
 */

#ifndef __GEOMETRY_TESSELATE_H
#define __GEOMETRY_TESSELATE_H

#include "../CGLA/Vec3f.h"
#include "../Util/Grid2D.h"
#include "IndexedFaceSet.h"

namespace Geometry
{
	
	class ParSurf
{
public:
	virtual CGLA::Vec3f operator()(float u, float v) const = 0;
};

extern float MAX_ERR;
extern float MAX_DIST;
extern int ADAPTIVE;


void tessellate(IndexedFaceSet& face_set, ParSurf& s,
			   float u_min, float u_max, float v_min, float v_max, 
			   int n, int m);

void tessellate(IndexedFaceSet& face_set, ParSurf& s, 
			   Util::Grid2D<CGLA::Vec3f>& inigrid);

void tessellate(IndexedFaceSet& face_set, ParSurf& s, 
			   std::vector<CGLA::Vec2f> uv_points,
			   std::vector<CGLA::Vec3i> triangles);
}

#endif
