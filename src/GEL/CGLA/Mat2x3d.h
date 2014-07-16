/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat2x3d.h
 * @brief 2x3 double matrix class
 */

#ifndef __CGLA_MAT2X3D_H__
#define __CGLA_MAT2X3D_H__

#include "Vec2d.h"
#include "Vec3d.h"
#include "ArithMatFloat.h"

namespace CGLA
{

	/**  \brief 2x3 double matrix class.

			 This class is useful for projecting a vector from 3D space to 2D.
	*/
	class Mat2x3d: public ArithMatFloat<Vec2d, Vec3d, Mat2x3d, 2>
	{

	public:
		/// Construct Mat2x3d from two Vec3f vectors (vectors become rows)
		Mat2x3d(const Vec3d& _a, const Vec3d& _b): 
			ArithMatFloat<Vec2d, Vec3d, Mat2x3d, 2> (_a,_b) {}

		/// Construct 0 matrix.
		Mat2x3d() {}
	};

	/**  \brief 3x2 double matrix class.

			 This class is useful for going from plane to 3D coordinates.
	*/
	class Mat3x2d: public ArithMatFloat<Vec3d, Vec2d, Mat3x2d, 3>
	{

	public:

		/** Construct matrix from three Vec2d vectors which become the 
				rows of the matrix. */
		Mat3x2d(const Vec2d& _a, const Vec2d& _b, const Vec2d& _c): 
			ArithMatFloat<Vec3d, Vec2d, Mat3x2d, 3> (_a,_b,_c) {}

		/// Construct 0 matrix.
		Mat3x2d() {}

	};


}
#endif
