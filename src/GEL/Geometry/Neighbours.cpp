/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Neighbours.h"

namespace Geometry
{
	using namespace CGLA;

	const Vec3f N6f[6] =
	{
		Vec3f(-1,0,0),
		Vec3f( 1,0,0),
		Vec3f( 0,-1,0),
		Vec3f( 0, 1,0),
		Vec3f( 0,0,-1),
		Vec3f( 0,0,1)
	};
	
	const Vec3i N6i[6] =
	{
		Vec3i(-1,0,0),
		Vec3i( 1,0,0),
		Vec3i( 0,-1,0),
		Vec3i( 0, 1,0),
		Vec3i( 0,0,-1),
		Vec3i( 0,0,1)		
	};

	const Vec3d N6d[6] =
	{
		Vec3d(-1,0,0),
		Vec3d( 1,0,0),
		Vec3d( 0,-1,0),
		Vec3d( 0, 1,0),
		Vec3d( 0,0,-1),
		Vec3d( 0,0,1)		
	};
	
	const Vec3f N26f[26] =
	{
		Vec3f(-1,-1,-1),
		Vec3f( 0,-1,-1),
		Vec3f( 1,-1,-1),
		Vec3f(-1, 0,-1),
		Vec3f( 0, 0,-1),
		Vec3f( 1, 0,-1),
		Vec3f(-1, 1,-1),
		Vec3f( 0, 1,-1),
		Vec3f( 1, 1,-1),

		Vec3f(-1,-1, 0),
		Vec3f( 0,-1, 0),
		Vec3f( 1,-1, 0),
		Vec3f(-1, 0, 0),
		Vec3f( 1, 0, 0),
		Vec3f(-1, 1, 0),
		Vec3f( 0, 1, 0),
		Vec3f( 1, 1, 0),

		Vec3f(-1,-1, 1),
		Vec3f( 0,-1, 1),
		Vec3f( 1,-1, 1),
		Vec3f(-1, 0, 1),
		Vec3f( 0, 0, 1),
		Vec3f( 1, 0, 1),
		Vec3f(-1, 1, 1),
		Vec3f( 0, 1, 1),
		Vec3f( 1, 1, 1)
		};

	const Vec3i N26i[26] = 
	{
		Vec3i(-1,-1,-1),
		Vec3i( 0,-1,-1),
		Vec3i( 1,-1,-1),
		Vec3i(-1, 0,-1),
		Vec3i( 0, 0,-1),
		Vec3i( 1, 0,-1),
		Vec3i(-1, 1,-1),
		Vec3i( 0, 1,-1),
		Vec3i( 1, 1,-1),

		Vec3i(-1,-1, 0),
		Vec3i( 0,-1, 0),
		Vec3i( 1,-1, 0),
		Vec3i(-1, 0, 0),
		Vec3i( 1, 0, 0),
		Vec3i(-1, 1, 0),
		Vec3i( 0, 1, 0),
		Vec3i( 1, 1, 0),

		Vec3i(-1,-1, 1),
		Vec3i( 0,-1, 1),
		Vec3i( 1,-1, 1),
		Vec3i(-1, 0, 1),
		Vec3i( 0, 0, 1),
		Vec3i( 1, 0, 1),
		Vec3i(-1, 1, 1),
		Vec3i( 0, 1, 1),
		Vec3i( 1, 1, 1)
	};

	const Vec3d N26d[26] = 
	{
		Vec3d(-1,-1,-1),
		Vec3d( 0,-1,-1),
		Vec3d( 1,-1,-1),
		Vec3d(-1, 0,-1),
		Vec3d( 0, 0,-1),
		Vec3d( 1, 0,-1),
		Vec3d(-1, 1,-1),
		Vec3d( 0, 1,-1),
		Vec3d( 1, 1,-1),

		Vec3d(-1,-1, 0),
		Vec3d( 0,-1, 0),
		Vec3d( 1,-1, 0),
		Vec3d(-1, 0, 0),
		Vec3d( 1, 0, 0),
		Vec3d(-1, 1, 0),
		Vec3d( 0, 1, 0),
		Vec3d( 1, 1, 0),

		Vec3d(-1,-1, 1),
		Vec3d( 0,-1, 1),
		Vec3d( 1,-1, 1),
		Vec3d(-1, 0, 1),
		Vec3d( 0, 0, 1),
		Vec3d( 1, 0, 1),
		Vec3d(-1, 1, 1),
		Vec3d( 0, 1, 1),
		Vec3d( 1, 1, 1)
	};


	const Vec3i CubeCorners8i[8] = 
		{
			Vec3i(0,0,0),
			Vec3i(1,0,0),
			Vec3i(0,1,0),
			Vec3i(1,1,0),
			Vec3i(0,0,1),
			Vec3i(1,0,1),
			Vec3i(0,1,1),
			Vec3i(1,1,1),
		};

	const Vec3f CubeCorners8f[8] = 
	{
		Vec3f(0,0,0),
		Vec3f(1,0,0),
		Vec3f(0,1,0),
		Vec3f(1,1,0),
		Vec3f(0,0,1),
		Vec3f(1,0,1),
		Vec3f(0,1,1),
		Vec3f(1,1,1),
	};
}
