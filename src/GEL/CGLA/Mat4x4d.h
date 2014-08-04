/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat4x4d.h
 * @brief 4x4 double matrix class
 */

#ifndef __CGLA_MAT4X4D_H__
#define __CGLA_MAT4X4D_H__

#include "ExceptionStandard.h"
#include "CGLA.h"
#include "Vec3d.h"
#include "Vec4d.h"
#include "ArithSqMat4x4Float.h"


namespace CGLA {


  /** \brief 4x4 double matrix.

      This class is useful for transformations such as perspective projections 
      or translation where 3x3 matrices do not suffice. */
  class Mat4x4d: public ArithSqMat4x4Float<Vec4d, Mat4x4d>
    {
    public:
  
      /// Construct a Mat4x4d from four Vec4d vectors
      Mat4x4d(Vec4d _a, Vec4d _b, Vec4d _c, Vec4d _d): 
	ArithSqMat4x4Float<Vec4d, Mat4x4d> (_a,_b,_c,_d) {}
  
      /// Construct the nan matrix
      Mat4x4d() {}

      /// Construct a matrix with identical elements.
      explicit Mat4x4d(double a): ArithSqMat4x4Float<Vec4d, Mat4x4d> (a) {}
    };

  /// Create a rotation _matrix. Rotates about one of the major axes.
  Mat4x4d rotation_Mat4x4d(CGLA::Axis axis, float angle);

  /// Create a translation matrix
  Mat4x4d translation_Mat4x4d(const Vec3d&);

  /// Create a scaling matrix.
  Mat4x4d scaling_Mat4x4d(const Vec3d&);

  /// Create an identity matrix.
  inline Mat4x4d identity_Mat4x4d()
    {
      return Mat4x4d(Vec4d(1.0,0.0,0.0,0.0), 
		     Vec4d(0.0,1.0,0.0,0.0), 
		     Vec4d(0.0,0.0,1.0,0.0), 
		     Vec4d(0.0,0.0,0.0,1.0));
    }

  /** Compute inverse assuming that the upper-left 3x3 sub-matrix is
      orthonormal (which is the case if the transformation is only
      a concatenation of rotations and translations).
  */
  inline Mat4x4d invert_ortho(const Mat4x4d& m)
  {
    Vec3d rx(m[0][0], m[1][0], m[2][0]);
    Vec3d ry(m[0][1], m[1][1], m[2][1]);
    Vec3d rz(m[0][2], m[1][2], m[2][2]);
    Vec3d t(m[0][3], m[1][3], m[2][3]);

    return Mat4x4d(Vec4d(rx, -dot(t, rx)),
		   Vec4d(ry, -dot(t, ry)),
		   Vec4d(rz, -dot(t, rz)),
		   Vec4d(0.0, 0.0, 0.0, 1.0));
  }   
}
#endif







