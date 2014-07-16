/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat3x3f.h
 * @brief 3x3 float matrix class
 */

#ifndef __CGLA_MAT3X3F_H__
#define __CGLA_MAT3X3F_H__

#include "CGLA.h"
#include "Vec3f.h"
#include "ArithSqMat3x3Float.h"

namespace CGLA 
{

  /** \brief 3 by 3 float matrix.

      This class will typically be used for rotation or
      scaling matrices for 3D vectors. */
  class Mat3x3f: public ArithSqMat3x3Float<Vec3f, Mat3x3f>
    {
    public:

      /// Construct matrix from 3 Vec3f vectors.
      Mat3x3f(Vec3f _a, Vec3f _b, Vec3f _c): 
	ArithSqMat3x3Float<Vec3f, Mat3x3f> (_a,_b,_c) {}
  
      /// Construct the 0 matrix
      Mat3x3f() {}

      /// Construct a matrix from a single scalar value.
      explicit Mat3x3f(float a): ArithSqMat3x3Float<Vec3f, Mat3x3f>(a) {}

    };

  /// Create a rotation _matrix. Rotates about one of the major axes.
  Mat3x3f rotation_Mat3x3f(CGLA::Axis axis, float angle);

  /// Create a scaling matrix.
  Mat3x3f scaling_Mat3x3f(const Vec3f&);

  /// Create an identity matrix.
  inline Mat3x3f identity_Mat3x3f()
    {
      return Mat3x3f(Vec3f(1.0f,0.0f,0.0f), Vec3f(0.0f,1.0f,0.0f), Vec3f(0.0f,0.0f,1.0f));
    }

}
#endif







