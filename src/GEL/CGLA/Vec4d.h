/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec4d.h
 * @brief 4D double vector class.
 */

#ifndef __CGLA_VEC4D_H__
#define __CGLA_VEC4D_H__

#include "Vec3d.h"
#include "ArithVec4Float.h"

namespace CGLA {

	/** \brief A four dimensional floating point vector. 

			This class is also used (via typedef) for
			homogeneous vectors.
	*/

	class Vec4d: public ArithVec4Float<double,Vec4d>
	{
	public:
  
		/// Construct a (0,0,0,0) homogenous Vector
		Vec4d(): ArithVec4Float<double,Vec4d>(0.0,0.0,0.0,0.0) {}

		/// Construct a (0,0,0,0) homogenous Vector
		explicit Vec4d(double _a): ArithVec4Float<double,Vec4d>(_a,_a,_a,_a) {}

		/// Construct a 4D vector
		Vec4d(double _a, double _b, double _c, double _d): 
			ArithVec4Float<double,Vec4d>(_a,_b,_c,_d) {}

		/// Construct a homogenous vector from a non-homogenous.
		explicit Vec4d(const Vec3d& v,double _d): 
			ArithVec4Float<double,Vec4d>(v[0],v[1],v[2],_d) {}
	};
}
#endif

