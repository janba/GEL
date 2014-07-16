/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "../CGLA/Vec3d.h"
#include "ThreeDDDA.h"

using namespace std;
using namespace CGLA;

namespace Geometry
{


	/** Return the position on the fine grid of a vector v. 
			We simply multiply by the precision and round down. Note that 
			it is necessary to use floor since simply clamping the value
			will always round toward zero.
	*/
	const Vec3i ThreeDDDA::grid_pos(const Vec3f& v) const
	{
		Vec3i g;
		for(int i=0;i<3;++i)
			g[i] = static_cast<int>(floor(v[i]*PRECISION));
		return g;
	}


	/// Compute the cell containing a given fine grid position.
	Vec3i ThreeDDDA::get_cell(const Vec3i& v) const
	{
		Vec3i c;
		for(int i=0;i<3;++i)
			{
				// If we are on the right side of the axis
				if(v[i]>=0)
					// We simply divide by the precision.
					c[i] = (v[i])/PRECISION;
				else
					// If we are on the left side, we have to add one
					// before and subtract one after dividing. This reflects
					// the lack of symmetry.
					c[i] = (v[i]+1)/PRECISION-1;
			}
		return c;
	}

	/// Mirror a fine grid position
	Vec3i ThreeDDDA::mirror(const Vec3i& v) const
	{
		Vec3i m;
		for(int i=0;i<3;++i)
			{
				// Mirroring simply inverts the sign, however, 
				// to reflect the fact that the cells on either side of the
				// axis are numbered differently, we have to add one before flipping 
				// the sign.
				if(sgn[i]==-1)
					m[i] = -(v[i]+1);
				else
					m[i] = v[i];
			}
		return m;
	}


	/** Create a 3DDDA based on the two end-points of the ray. An optional
			last argument indicates the resolution of the grid. */
	ThreeDDDA::ThreeDDDA(const CGLA::Vec3f& p0, const CGLA::Vec3f& p1, 
											 int _PRECISION):
		PRECISION(_PRECISION),
		sgn(p1[0]>=p0[0] ? 1 : -1,
				p1[1]>=p0[1] ? 1 : -1,
				p1[2]>=p0[2] ? 1 : -1),
		p0_int(grid_pos(p0)),
		p1_int(grid_pos(p1)),
		dir(sgn*(p1_int-p0_int)),
		first_cell(get_cell(p0_int)),
		last_cell(get_cell(p1_int)),
		no_steps(dot(sgn, last_cell-first_cell)),
		dir_pre_mul(2*dir*PRECISION)
	{
		// Compute the mirrored position of the ray starting point
		const Vec3i p0_int_m = mirror(p0_int);

		// Compute the cell of the mirrored position.
		const Vec3i first_cell_m = get_cell(p0_int_m);

		/* Finally, compute delta = p - p0 where 
			 "p" is the position of the far,upper,right corner of the cell containing
			 the ray starting point and "p0" is the starting point (mirrored). 
			 Since the ray starts at a grid position + (0.5,0.5,0.5) we have to
			 add one half to each grid coordinate of p0. Since we are doing integer
			 arithmetic, we simply multiply everything by two instead.
		*/
		const Vec3i delta = 2*(first_cell_m*PRECISION+Vec3i(PRECISION))-
			(2*p0_int_m  + Vec3i(1));

		/* Compute the discriminators. These are (x-x0) * dir_y, (y-y0)*dir_x
			 and so on for all permutations. These will be updated in the course
			 of the algorithm. */ 
		for(int i=0;i<3;++i)
			{
				const int k1 = (i + 1) % 3;
				const int k2 = (i + 2) % 3;
				disc[i][k1] = static_cast<Integer64Bits>(delta[i])*dir[k1];
				disc[i][k2] = static_cast<Integer64Bits>(delta[i])*dir[k2];
			}

		// Finally, we set the cell position equal to the first_cell and we
		// are ready to iterate.
		cell = first_cell;
	}

	/// Get the initial point containing the ray.
	const CGLA::Vec3f ThreeDDDA::step() 
	{
		int step_dir;
		if(disc[1][0] > disc[0][1])
			if(disc[2][0] > disc[0][2])
				step_dir = 0;
			else
				step_dir = 2;
		else
			if(disc[2][1] > disc[1][2])
				step_dir = 1;
			else
				step_dir = 2;
		
		const int k1 = (step_dir + 1) % 3;
		const int k2 = (step_dir + 2) % 3;
		
		double delta;
		delta=cell[step_dir]+sgn[step_dir]-double(p0_int[step_dir]+0.5f)/PRECISION;
		Vec3f real_dir = get_last_point() - get_first_point();
		double ck1 = real_dir[k1]/real_dir[step_dir];
		double ck2 = real_dir[k2]/real_dir[step_dir];

		Vec3f p;
		p[step_dir] = delta;
		p[k1] = ck1 * delta;
		p[k2] = ck2 * delta;
 		p += get_first_point();

		cell[step_dir] += sgn[step_dir];
		disc[step_dir][k1] += dir_pre_mul[k1];
		disc[step_dir][k2] += dir_pre_mul[k2];
		return Vec3f(p);
	}

	
}
