/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat2x2f.h
 * @brief 2x2 float matrix class
 */

#ifndef __CGLA_MAT2X2F_H__
#define __CGLA_MAT2X2F_H__

#include "Vec2f.h"
#include "ArithSqMat2x2Float.h"


namespace CGLA 
{

  /** \brief Two by two float matrix. 

	This class is useful for various 
	vector transformations in the plane. */
  class Mat2x2f: public ArithSqMat2x2Float<Vec2f, Mat2x2f>
    {
    public:

      /// Construct a Mat2x2f from two Vec2f vectors.
      Mat2x2f(Vec2f _a, Vec2f _b): ArithSqMat2x2Float<Vec2f, Mat2x2f> (_a,_b) {}

      /// Construct a Mat2x2f from four scalars.
      Mat2x2f(float _a, float _b, float _c, float _d): 
	ArithSqMat2x2Float<Vec2f, Mat2x2f>(Vec2f(_a,_b),Vec2f(_c,_d)) {}
  
      /// Construct the NAN matrix
      Mat2x2f() {}
      
      /// Construct a Mat2x2f from a single scalar
      explicit Mat2x2f(float a): ArithSqMat2x2Float<Vec2f, Mat2x2f>(a) {}


    };

}
#endif
