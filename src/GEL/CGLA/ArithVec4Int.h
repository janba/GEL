/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract 4D integer vector class
 * ----------------------------------------------------------------------- */

/** @file ArithVec4Int.h
 * @brief Abstract 4D integer vector class
 */

#ifndef CGLA__ARITHVEC4INT_H
#define CGLA__ARITHVEC4INT_H

#include <GEL/CGLA/ArithVecInt.h>

namespace CGLA {

	template<class T, class V>
	class ArithVec4Int: public ArithVecInt<T,V,4>
	{
	public:

		/// Construct a 4D int vector.
		constexpr ArithVec4Int(T a, T b, T c, T d): ArithVecInt<T,V,4>(a,b,c,d) {}

		/// Construct a 4D int vector.
		constexpr ArithVec4Int() = default;
		
	};
}

#endif

