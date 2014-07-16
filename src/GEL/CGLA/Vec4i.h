/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec4i.h
 * @brief 4D integer vector class.
 */

#ifndef __CGLA_VEC4I_H__
#define __CGLA_VEC4I_H__

#include "ArithVec4Int.h"

namespace CGLA 
{
	class Vec4f;
	class Vec4uc;
	class Vec4usi;

	/** \brief 4D integer vector. 

	    This class does not really extend the template
			and hence provides only the basic facilities of an ArithVec. 
			The class is typically used for indices to 4D voxel grids. */
	class Vec4i: public ArithVec4Int<int,Vec4i>
	{
	public:
  
		/// Construct 0 vector.
		Vec4i() {}

		/// Construct a 4D integer vector.
		Vec4i(int _a,int _b,int _c, int _d): ArithVec4Int<int,Vec4i>(_a,_b,_c,_d) {}

		/// Construct a 4D integer vector with 4 identical coordinates.
		explicit Vec4i(int a): ArithVec4Int<int,Vec4i>(a,a,a,a) {}
	
		/// Construct from a Vec4f.
		explicit Vec4i(const Vec4f& v);

		/// Construct from a Vec4uc.
		explicit Vec4i(const Vec4uc& v);

		/// Construct from a Vec4usi.
		explicit Vec4i(const Vec4usi& v);

	};


}
#endif
