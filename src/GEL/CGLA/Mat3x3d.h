/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat3x3d.h
 * @brief 3x3 double matrix class
 */

#ifndef __CGLA_MAT3X3D_H__
#define __CGLA_MAT3X3D_H__

#include "CGLA.h"
#include "Vec3d.h"
#include "ArithSqMat3x3Float.h"

namespace CGLA {

	/** \brief 3 by 3 double matrix.

			This class will typically be used for rotation or
			scaling matrices for 3D vectors. */
	class Mat3x3d: public ArithSqMat3x3Float<Vec3d, Mat3x3d>
	{
	public:

		/// Construct matrix from 3 Vec3d vectors.
		Mat3x3d(Vec3d _a, Vec3d _b, Vec3d _c): 
			ArithSqMat3x3Float<Vec3d, Mat3x3d> (_a,_b,_c) {}
  
		/// Construct the 0 matrix
		Mat3x3d() {}

		/// Construct a matrix from a single scalar value.
		explicit Mat3x3d(float a): ArithSqMat3x3Float<Vec3d, Mat3x3d>(a) {}

	};

  /// Create a rotation _matrix. Rotates about one of the major axes.
  Mat3x3d rotation_Mat3x3d(CGLA::Axis axis, double angle);

  /// Create a scaling matrix.
  Mat3x3d scaling_Mat3x3d(const Vec3d&);


	/// Create an identity matrix.
	inline Mat3x3d identity_Mat3x3d()
	{
		return Mat3x3d(Vec3d(1.0,0.0,0.0), Vec3d(0.0,1.0,0.0), Vec3d(0.0,0.0,1.0));
	}

}
#endif







