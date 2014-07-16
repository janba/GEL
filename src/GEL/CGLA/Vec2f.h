/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec2f.h
 * @brief 2D float vector class.
 */

#ifndef __CGLA_VEC2F_H__
#define __CGLA_VEC2F_H__

#include "ArithVec2Float.h"
#include "Vec2i.h"


namespace CGLA {

	/** \brief 2D floating point vector */

	class Vec2f: public ArithVec2Float<float,Vec2f>
	{
	public:

		Vec2f() {}
		Vec2f(float _a,float _b): ArithVec2Float<float,Vec2f>(_a,_b) {}

		template<class T, class V, unsigned int N>
		explicit Vec2f(const ArithVec<T,V,N>& v): 
			ArithVec2Float<float,Vec2f>(static_cast<float>(v[0]),
																	static_cast<float>(v[1])) {}

		explicit Vec2f(float a): ArithVec2Float<float,Vec2f>(a,a) {}
	};



}
#endif
