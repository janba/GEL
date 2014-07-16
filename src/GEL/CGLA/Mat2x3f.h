/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat2x3f.h
 * @brief 2x3 float matrix class
 */

#ifndef __CGLA_MAT2X3F_H__
#define __CGLA_MAT2X3F_H__

#include "Vec2f.h"
#include "Vec3f.h"
#include "ArithMatFloat.h"

namespace CGLA
{

  /**  \brief 2x3 float matrix class.

       This class is useful for projecting a vector from 3D space to 2D.
  */
  class Mat2x3f: public ArithMatFloat<Vec2f, Vec3f, Mat2x3f, 2>
    {

    public:
      /// Construct Mat2x3f from two Vec3f vectors (vectors become rows)
      Mat2x3f(const Vec3f& _a, const Vec3f& _b): 
	ArithMatFloat<Vec2f, Vec3f, Mat2x3f, 2> (_a,_b) {}

      /// Construct NAN matrix.
      Mat2x3f() {}
    };


  /**  \brief 3x2 float matrix class.

       This class is useful for going from plane to 3D coordinates.
  */
  class Mat3x2f: public ArithMatFloat<Vec3f, Vec2f, Mat3x2f, 3>
    {

    public:

      /** Construct matrix from three Vec2f vectors which become the 
	  rows of the matrix. */
      Mat3x2f(const Vec2f& _a, const Vec2f& _b, const Vec2f& _c): 
	ArithMatFloat<Vec3f, Vec2f, Mat3x2f, 3> (_a,_b,_c) {}

      /// Construct NAN matrix.
      Mat3x2f() {}
    };

}
#endif
