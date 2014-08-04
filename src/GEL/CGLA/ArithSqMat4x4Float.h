/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract 4x4 floating point matrix class
 * ----------------------------------------------------------------------- */

/** @file ArithSqMat4x4Float.h
 * @brief Abstract 4x4 floating point matrix class
 */

#ifndef __CGLA_ARITHSQMAT4X4FLOAT_H
#define __CGLA_ARITHSQMAT4X4FLOAT_H

#include "ExceptionStandard.h"
#include "CGLA.h"
#include "Vec3f.h"
#include "Vec4f.h"
#include "ArithSqMatFloat.h"


namespace CGLA 
{
  CGLA_DERIVEEXCEPTION(Mat4x4fException)
  CGLA_DERIVEEXCEPTION(Mat4x4fNotAffine)
  CGLA_DERIVEEXCEPTION(Mat4x4fSingular)
    
  /** \brief 4 by 4 float matrix template.
			
      this class template is useful for transformations such as perspective 
      projections or translation where 3x3 matrices do not suffice. */
  template<class VT, class M>
  class ArithSqMat4x4Float: public ArithSqMatFloat<VT, M, 4>
  {
  public:
    
    /// Vector type
    typedef VT VectorType;
    
    /// The type of a matrix element
    typedef typename VT::ScalarType ScalarType;
    
  public:
    
    /// Construct a Mat4x4f from four V vectors
    ArithSqMat4x4Float(VT a, VT b, VT c, VT d): 
      ArithSqMatFloat<VT, M, 4> (a,b,c,d) {}
  
    /// Construct the NAN matrix
    ArithSqMat4x4Float() {}

    /// Construct matrix where all values are equal to constructor argument.
    explicit ArithSqMat4x4Float(ScalarType  _a):
      ArithSqMatFloat<VT,M,4>(_a) {}

    /** Multiply vector onto matrix. Here the fourth coordinate 
	is se to 0. This removes any translation from the matrix.
	Useful if one wants to transform a vector which does not
	represent a point but a direction. Note that this is not
	correct for transforming normal vectors if the matric 
	contains anisotropic scaling. */
    template<class T, class VecT>
    const VecT mul_3D_vector(const ArithVec3Float<T,VecT>& v_in) const
    {
      VT v_out  = (*this) * VT(v_in[0],v_in[1],v_in[2],0);
      return VecT(v_out[0],v_out[1],v_out[2]);
    }

    /** Multiply 3D point onto matrix. Here the fourth coordinate 
	becomes 1 to ensure that the point is translated. Note that
	the vector is converted back into a Vec3f without any 
	division by w. This is deliberate: Typically, w=1 except
	for projections. If we are doing projection, we can use
	project_3D_point instead */
    template<class T, class VecT>
    const VecT mul_3D_point(const ArithVec3Float<T,VecT> & v_in) const
    {
      VT v_out  = (*this) * VT(v_in[0],v_in[1],v_in[2],1);
      return VecT(v_out[0],v_out[1],v_out[2]);
    }

    /** Multiply 3D point onto matrix. We set w=1 before 
	multiplication and divide by w after multiplication. */
    template<class T, class VecT>
    const VecT project_3D_point(const ArithVec3Float<T,VecT>& v_in) const
    {
      VT v_out = (*this) * VT(v_in[0],v_in[1],v_in[2],1);
      v_out.de_homogenize();
      return VecT(v_out);
    }

  };

  /** Compute the adjoint of a matrix. This is the matrix where each
      entry is the subdeterminant of 'in' where the row and column of 
      the element is removed. Use mostly to compute the inverse */
  template<class VT, class M>
  M adjoint(const ArithSqMat4x4Float<VT,M>&);
		
	/** Compute the determinant of a 4x4 matrix. The code below is what
			I found to be most robust. The original implementation used direct
			computation of the 3x3 sub-determinants and I also tried a direct
			computation based on enumerating all permutations. The code below
			is better. */
	template<class V, class M>
			inline double determinant(const ArithSqMat4x4Float<V,M>& m)
			{
/* 					return m[0][0] * (m[1][1]*(m[2][2]*m[3][3]-m[3][2]*m[2][3])- */
/* 														m[2][1]*(m[1][2]*m[3][3]-m[3][2]*m[1][3])+ */
/* 														m[3][1]*(m[1][2]*m[2][3]-m[2][2]*m[1][3])) */
/* 							- m[1][0] * (m[0][1]*(m[2][2]*m[3][3]-m[3][2]*m[2][3])- */
/* 													 m[2][1]*(m[0][2]*m[3][3]-m[3][2]*m[0][3])+ */
/* 													 m[3][1]*(m[0][2]*m[2][3]-m[2][2]*m[0][3])) */
/* 							+ m[2][0] * (m[0][1]*(m[1][2]*m[3][3]-m[3][2]*m[1][3])- */
/* 													 m[1][1]*(m[0][2]*m[3][3]-m[3][2]*m[0][3])+ */
/* 													 m[3][1]*(m[0][2]*m[1][3]-m[1][2]*m[0][3])) */
/* 							- m[3][0] * (m[0][1]*(m[1][2]*m[2][3]-m[2][2]*m[1][3])- */
/* 													 m[1][1]*(m[0][2]*m[2][3]-m[2][2]*m[0][3])+ */
/* 													 m[2][1]*(m[0][2]*m[1][3]-m[1][2]*m[0][3])); */
		typedef typename M::ScalarType ScalarType;
    ScalarType a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;
    
    /* assign to individual variable names to aid selecting */
		/*  correct elements */
		
		a1 = m[0][0]; b1 = m[1][0]; c1 = m[2][0]; d1 = m[3][0];
		a2 = m[0][1]; b2 = m[1][1]; c2 = m[2][1]; d2 = m[3][1];
		a3 = m[0][2]; b3 = m[1][2]; c3 = m[2][2]; d3 = m[3][2];
		a4 = m[0][3]; b4 = m[1][3]; c4 = m[2][3]; d4 = m[3][3];

					return 
							a1 * (b2*(c3*d4-d3*c4)-c2*(b3*d4-d3*b4)+d2*(b3*c4-c3*b4))
							- b1 * (a2*(c3*d4-d3*c4)-c2*(a3*d4-d3*a4)+d2*(a3*c4-c3*a4))
							+ c1 * (a2*(b3*d4-d3*b4)-b2*(a3*d4-d3*a4)+d2*(a3*b4-b3*a4))
							- d1 * (a2*(b3*c4-c3*b4)-b2*(a3*c4-c3*a4)+c2*(a3*b4-b3*a4));

			}

  /// Compute the inverse matrix of a Mat4x4f.
  template<class VT, class M>
  M invert(const ArithSqMat4x4Float<VT,M>&);

  /// Compute the inverse matrix of a Mat4x4f that is affine
  template<class VT, class M>
  M invert_affine(const ArithSqMat4x4Float<VT,M>&);

}
#endif







