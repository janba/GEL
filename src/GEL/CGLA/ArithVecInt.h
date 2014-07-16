/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithVecInt.h
 * @brief Abstract integer vector class
 */

#ifndef __CGLA__ARITHVECINT_H__
#define __CGLA__ARITHVECINT_H__

#include "ArithVec.h"

namespace CGLA {

	template<class T, class V, unsigned int N>
	class ArithVecInt: public ArithVec<T,V,N>
	{
	public:

		ArithVecInt() {}

		ArithVecInt(T a): 
			ArithVec<T,V,N>(a) {}

		ArithVecInt(T a, T b): 
			ArithVec<T,V,N>(a,b) {}

		ArithVecInt(T a, T b, T c): 
			ArithVec<T,V,N>(a,b,c) {}

		ArithVecInt(T a, T b, T c, T d): 
			ArithVec<T,V,N>(a,b,c,d) {}

	};
}

#endif

