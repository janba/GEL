/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithVec3Int.h
  * Abstract 3D integer vector class.
  */

#ifndef __CGLA__ARITHVEC3INT_H__
#define __CGLA__ARITHVEC3INT_H__

#include "ArithVecInt.h"

namespace CGLA {

	template<class T, class V>
	class ArithVec3Int: public ArithVecInt<T,V,3>
	{
	public:

		/// Construct a 3D int vector.
		ArithVec3Int(T a, T b, T c): ArithVecInt<T,V,3>(a,b,c) {}

		/// Construct a 3D int vector.
		ArithVec3Int() {}
		
	};

	/// Returns cross product of arguments
	template<class T, class V>
	inline V cross( const ArithVec3Int<T,V>& x, 
									const ArithVec3Int<T,V>& y ) 
	{
		return V( x[1] * y[2] - x[2] * y[1], 
							x[2] * y[0] - x[0] * y[2], 
							x[0] * y[1] - x[1] * y[0] );
	}

}

#endif

