/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithSqMat3x3Float.h
 * @brief Abstract 3x3 floating point matrix class
 */

#ifndef __CGLA_ARITHSQMAT3X3FLOAT_H__
#define __CGLA_ARITHSQMAT3X3FLOAT_H__

#include "ExceptionStandard.h"
#include "CGLA.h"
#include "Vec3f.h"
#include "ArithSqMatFloat.h"

namespace CGLA 
{

  CGLA_DERIVEEXCEPTION(Mat3x3fException)
    CGLA_DERIVEEXCEPTION(Mat3x3fSingular)

    /** \brief 3 by 3 float matrix template.
				
		This class template will typically be used for rotation or
		scaling matrices for 3D vectors. */
    template<class V, class M>
    class ArithSqMat3x3Float: public ArithSqMatFloat<V, M, 3>
    {
    public:

      /// Vector type
      typedef V VectorType;

      /// The type of a matrix element
      typedef typename V::ScalarType ScalarType;

    public:

      /// Construct matrix from 3 Vec3f vectors.
      ArithSqMat3x3Float(V _a, V _b, V _c): 
	ArithSqMatFloat<V, M, 3> (_a,_b,_c) {}
  
      /// Construct the 0 matrix
      ArithSqMat3x3Float() {}

      /// Construct a matrix from a single scalar value.
      explicit ArithSqMat3x3Float(ScalarType a): 
	ArithSqMatFloat<V, M, 3>(a) {}

    };

  /// Invert 3x3 matrix
  template<class V, class M>
    M invert(const ArithSqMat3x3Float<V,M>&);

  /** Compute determinant. There is a more generic function for
      computing determinants of square matrices (ArithSqMat). This one
      is faster but works only on Mat3x3f */
  template<class V, class M>
    inline 
    typename ArithSqMat3x3Float<V,M>::ScalarType
    determinant(const ArithSqMat3x3Float<V,M>& m)
    {
      return 
	m[0][0]*m[1][1]*m[2][2] +
	m[0][1]*m[1][2]*m[2][0] +
	m[0][2]*m[1][0]*m[2][1] -
	m[0][2]*m[1][1]*m[2][0] -
	m[0][0]*m[1][2]*m[2][1] -
	m[0][1]*m[1][0]*m[2][2] ;
    }


}
#endif







