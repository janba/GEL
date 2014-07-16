/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Vec3i.h
 * 3D integer vector class.
 */

#ifndef __CGLA_VEC3I_H__
#define __CGLA_VEC3I_H__

#include "ArithVec3Int.h"

namespace CGLA 
{
	class Vec3f;
	class Vec3d;
	class Vec3uc;
	class Vec3usi;

	/** \brief 3D integer vector. 

	    This class does not really extend the template
			and hence provides only the basic facilities of an ArithVec. 
			The class is typically used for indices to 3D voxel grids. */
	class Vec3i: public ArithVec3Int<int,Vec3i>
	{
	public:
  
		/// Construct 0 vector.
		Vec3i() {}

		/// Construct a 3D integer vector.
		Vec3i(int _a,int _b,int _c): ArithVec3Int<int,Vec3i>(_a,_b,_c) {}

		/// Construct a 3D integer vector with 3 identical coordinates.
		explicit Vec3i(int a): ArithVec3Int<int,Vec3i>(a,a,a) {}
	
		/// Construct from a Vec3f.
		explicit Vec3i(const Vec3f& v);

		/// Construct from a Vec3f.
		explicit Vec3i(const Vec3d& v);

		/// Construct from a Vec3uc.
		explicit Vec3i(const Vec3uc& v);

		/// Construct from a Vec3usi.
		explicit Vec3i(const Vec3usi& v);

	};


}
#endif
