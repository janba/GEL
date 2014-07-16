/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3d.h
 * @brief 3D double vector class.
 */

#ifndef __CGLA_VEC3D_H__
#define __CGLA_VEC3D_H__

#include "ArithVec3Float.h"
#include "Vec3i.h"
#include "Vec3usi.h"
#include "Vec3f.h"


namespace CGLA {

	/** \brief A 3D double vector. 

	Useful for high precision arithmetic. */

	class Vec3d: public ArithVec3Float<double,Vec3d>
	{
	public:

		/// Construct 0 vector
		Vec3d(){}

		/// Construct vector
		Vec3d(double a, double b, double c): ArithVec3Float<double,Vec3d>(a,b,c) {}

		/// Construct vector where all coords = a 
		explicit Vec3d(double a):
			ArithVec3Float<double,Vec3d>(a,a,a) {}

		/// Convert from int vector
		explicit Vec3d(const Vec3i& v): 
			ArithVec3Float<double,Vec3d>(v[0],v[1],v[2]) {}

		/// Construct from a 3D unsigned int vector.
		explicit Vec3d(const Vec3usi& v): 
			ArithVec3Float<double,Vec3d>(v[0],v[1],v[2]) {}

		/// Convert from float vector
		explicit Vec3d(const Vec3f& v): 
			ArithVec3Float<double,Vec3d>(v[0],v[1],v[2]) {}
	};


}
#endif
