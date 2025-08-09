/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Neighbours.h
 * @brief A number of arrays with relative 3D voxel indices.
 */

#ifndef GEOMETRY_NEIGHBOURS_H
#define GEOMETRY_NEIGHBOURS_H

#include <GEL/CGLA/Vec.h>

namespace Geometry
{
	extern const CGLA::Vec3f N6f[6];
	extern const CGLA::Vec3i N6i[6];
	extern const CGLA::Vec3d N6d[6];
	
	extern const CGLA::Vec3f N26f[26];
	extern const CGLA::Vec3i N26i[26];
	extern const CGLA::Vec3d N26d[26];

	extern const CGLA::Vec3i CubeCorners8i[8];
	extern const CGLA::Vec3f CubeCorners8f[8];
}


#endif
