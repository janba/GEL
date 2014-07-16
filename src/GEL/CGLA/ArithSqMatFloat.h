/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract floating point matrix class
 * ----------------------------------------------------------------------- */

/** @file ArithSqMatFloat.h
 * @brief Abstract floating point matrix class
 */

#ifndef __CGLA_ARITHSQMATFLOAT_H__
#define __CGLA_ARITHSQMATFLOAT_H__

#include "ArithMatFloat.h"

namespace CGLA 
{

  /** Template for square matrices.

      Some functions like trace and determinant work only on
      square matrices. To express this in the class hierarchy,
      ArithSqMatFloat was created. ArithSqMatFloat is derived from ArithMat
      and contains a few extra facilities applicable only to
      square matrices.
  */
  template <class VT, class MT, unsigned int ROWS>
    class ArithSqMatFloat: public ArithMatFloat<VT,VT,MT,ROWS> 
    { 
    public:

      /// Vector type
      typedef VT VectorType;

      /// The type of a matrix element
      typedef typename VT::ScalarType ScalarType;

    protected:

      /// Construct 0 matrix
      ArithSqMatFloat() {}

      /// Construct matrix where all values are equal to constructor argument.
      explicit ArithSqMatFloat(ScalarType _a):
	ArithMatFloat<VT,VT,MT,ROWS>(_a) {}

      /// Construct 2x2 Matrix from two vectors
      ArithSqMatFloat(VT _a, VT _b): 
	ArithMatFloat<VT,VT,MT,ROWS>(_a,_b) {}

      /// Construct 3x3 Matrix from three vectors
      ArithSqMatFloat(VT _a, VT _b, VT _c): 
	ArithMatFloat<VT,VT,MT,ROWS>(_a,_b,_c) {}

      /// Construct 4x4 Matrix from four vectors
      ArithSqMatFloat(VT _a, VT _b, VT _c, VT _d): 
	ArithMatFloat<VT,VT,MT,ROWS>(_a,_b,_c,_d) {}
		
    public:

      /** Assignment multiplication of matrices. 
	  This function is not very efficient. This because we need a temporary
	  matrix anyway, so it can't really be made efficient. */
      const MT& operator*=(const MT& m2)
	{
	  (*this) = (*this) * m2;
	  return static_cast<const MT&>(*this);
	}
		
      const MT& operator *=(ScalarType k) 
	{
	  return ArithMatFloat<VT,VT,MT,ROWS>::operator*=(k);
	}

      void identity()
	{
	  for(unsigned int i=0;i<ROWS;++i)
	    {
	      for(unsigned int j=0;j<ROWS;++j)
		(*this)[i][j] = ScalarType(0);
	      (*this)[i][i] = ScalarType(1);
	    }
	}

		
    };

  /** Multiply two matrices derived from same type, 
      producing a new of same type */
  template <class VT, class MT, unsigned int ROWS>
    inline MT operator*(const ArithSqMatFloat<VT,MT,ROWS>& m1,
			const ArithSqMatFloat<VT,MT,ROWS>& m2) 
    {
      MT n;
      for(unsigned int i=0;i<ROWS;i++)
	for(unsigned int j=0;j<ROWS;j++)
	  {
	    n[i][j] = 0;
	    for(unsigned int k=0;k<ROWS;k++)
	      n[i][j] += m1[i][k] * m2[k][j]; 
	  }
      return n;
    }

  /** Compute the transpose of a square matrix. This function returns
      the transpose of its argument. */
  template <class VT, class MT, unsigned int ROWS>
    inline MT transpose(const ArithSqMatFloat<VT,MT,ROWS>& m) 
    {
      MT m_new;
      for(unsigned int i=0;i<MT::get_v_dim();i++)
	for(unsigned int j=0;j<MT::get_h_dim();j++)
	  m_new[i][j] = m[j][i];
      return m_new;
    }

  /// Compute trace. Works only for sq. matrices.
  template <class VT, class MT, unsigned int ROWS>
    inline typename MT::ScalarType trace(const ArithSqMatFloat<VT,MT,ROWS>& M)
    {
      typename ArithSqMatFloat<VT,MT,ROWS>::ScalarType s=0;
      for(unsigned int i=0;i<ROWS;i++)
	s += M[i][i];
      return s;
    }

}
#endif
