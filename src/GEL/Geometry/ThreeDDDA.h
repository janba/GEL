/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file ThreeDDDA.h
 * @brief a class for traversing cells in a regular grid pierced by a ray.
 */
#ifndef __GEOMETRY_THREEDDDA_H
#define __GEOMETRY_THREEDDDA_H

#ifdef _MSC_VER 
typedef __int64 Integer64Bits;
#else
typedef int64_t Integer64Bits;
#endif


#include "../CGLA/Vec3f.h"
#include "../CGLA/Vec3i.h"

namespace Geometry
{


	/** \brief A ThreeDDDA is a class for traversing a grid of cells. 

      We wish to 
			enumerate all cells pierced by the ray going from a point p0 to a 
			point p1. It is dangerous to use floating point arithmetic since rounding
			errors may accumulate entailing that the 3ddda misses the cell containing
			the end point of the ray.
		
			Daniel Cohen-Or devised a technique for exact ray traversal using only 
			integer arithmetic. Using integer arithmetic, the computations are
			exact, but the ray must start and terminate exactly at grid points.
			I.e. the ray must start and terminate in cell corners.

			To overcome this issue, I have implemented a scheme where the ray
			starts and terminates on a grid that is finer than the cell grid.
			This means that we have fast and exact computations, and the ray can
			travel between two almost arbitrary points.

			A central problem with this scheme was that - in order to simplify
			things during the iterative traversal - the ray direction vector
			is always mirrored into the first octant (+x +y +z).  If the fine
			grid used for the ray start and end points is a superset of the
			grid points of the cell grid, this mirroring is almost impossible
			to implement correctly. This is due to the fact that the first
			cell on the right side of the y axis is called 0,y whereas the
			first cell on the left is called -1,y. This lack of symmetry makes
			things very difficult, but by ensuring that the endpoints of a ray
			can never lie on an a corner, edge, or face of the cell grid, we
			can make the problems go away.

			Hence, in this implementation, the fine grid is obtained by dividing 
			the cell into NxNxN subcells, but the grid points	of the fine grid lie 
			at the CENTRES of these subcells.  
		
	*/

	class ThreeDDDA
	{
		/** Resolution of the grid. This number indicates how many subcells
				(along each axis) a cell is divided into. */
		const int PRECISION;

		/// Sign indicates which octant contains the direction vector.
		const CGLA::Vec3i sgn;

		/// Fine grid position where the ray begins
		const CGLA::Vec3i p0_int;

		// Fine grid position where the ray ends.
		const CGLA::Vec3i p1_int;

		/// The direction of the ray 
		const CGLA::Vec3i dir;

		/// The cell containing the initial point on the ray.
		const CGLA::Vec3i first_cell;

		/// The cell containing the last point on the ray.
		const CGLA::Vec3i last_cell;

		/** Number of steps. We can compute the number of steps that the
				ray will take in advance. */
		const int no_steps; 

		/** The direction vector, premultiplied by the step length.
				This constant is used to update the equations governing the
				traversal. */
		const CGLA::Vec3i dir_pre_mul;

		/** Discriminators. These values are delta[i]*dir[j] where delta[i]
				is the distance from the origin of the ray to the current position
				along the i axis multiplied by two). dir[j] is the direction vector
				coordinate j. */
	
		Integer64Bits disc[3][3];

		/** The current cell containing the ray. This value along with the 
				discriminators mentioned above represent the state of a ThreeDDDA. */
		CGLA::Vec3i cell;


		/** Return the position on the fine grid of a vector v. 
				We simply multiply by the precision and round down. Note that 
				it is necessary to use floor since simply clamping the value
				will always round toward zero.
		*/
		const CGLA::Vec3i grid_pos(const CGLA::Vec3f& v) const;
	
		/// Compute the cell containing a given fine grid position.
		CGLA::Vec3i get_cell(const CGLA::Vec3i& v) const;
	
		/// Mirror a fine grid position
		CGLA::Vec3i mirror(const CGLA::Vec3i& v) const;


	public:

		/** Create a 3DDDA based on the two end-points of the ray. An optional
				last argument indicates the resolution of the grid. */
		ThreeDDDA(const CGLA::Vec3f& p0, const CGLA::Vec3f& p1, 
							int _PRECISION = 0x100);
	
		/// Return the total number of steps the ray will traverse.
		int get_no_steps() const
		{
			return no_steps;
		}

		/// Yields the current cell containing the ray.
		const CGLA::Vec3i& get_current_cell() const
		{
			return cell;
		}

		/// Get the initial cell containing the ray.
		const CGLA::Vec3i& get_first_cell() const
		{
			return first_cell;
		}

		/// Get the last cell containing the ray.
		const CGLA::Vec3i& get_last_cell() const
		{
			return last_cell;
		}

		/// Get the initial point containing the ray.
		const CGLA::Vec3f get_first_point() const
		{
			return (CGLA::Vec3f(p0_int)+CGLA::Vec3f(0.5))/PRECISION;
		}


		/// Get the last cell containing the ray.
		const CGLA::Vec3f get_last_point() const
		{
			return (CGLA::Vec3f(p1_int)+CGLA::Vec3f(0.5))/PRECISION;
		}
	
		/// Returns true if we have reached the last cell.
		bool at_end() const
		{
			return cell == last_cell;
		}

		/// Take a step along the ray.
		void operator++()
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
			cell[step_dir] += sgn[step_dir];
			disc[step_dir][k1] += dir_pre_mul[k1];
			disc[step_dir][k2] += dir_pre_mul[k2];
		}

		/// Get the initial point containing the ray.
			const CGLA::Vec3f step();


	};
}

#endif
