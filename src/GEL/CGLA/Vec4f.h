/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec4f.h
 * @brief 4D float vector class.
 */

#ifndef __CGLA_VEC4F_H__
#define __CGLA_VEC4F_H__

#include "Vec3f.h"
#include "ArithVec4Float.h"

namespace CGLA {

	/** \brief A four dimensional floating point vector. 

			This class is also used (via typedef) for
			homogeneous vectors.
	*/

	class Vec4f: public ArithVec4Float<float,Vec4f>
	{
	public:
  
		/// Construct a (0,0,0,0) homogenous Vector
		Vec4f(): ArithVec4Float<float,Vec4f>(0.0f,0.0f,0.0f,0.0f) {}

		/// Construct a (0,0,0,0) homogenous Vector
		explicit Vec4f(float _a): ArithVec4Float<float,Vec4f>(_a,_a,_a,_a) {}

		/// Construct a 4D vector
		Vec4f(float _a, float _b, float _c, float _d): 
			ArithVec4Float<float,Vec4f>(_a,_b,_c,_d) {}

		/// Construct a homogenous vector from a non-homogenous.
		explicit Vec4f(const Vec3f& v,float _d): 
			ArithVec4Float<float,Vec4f>(v[0],v[1],v[2],_d) {}


	};
}
#endif

