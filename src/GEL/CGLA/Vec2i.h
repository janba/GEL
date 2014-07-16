/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec2i.h
 * @brief 2D int vector class.
 */

#ifndef __CGLA_VEC2I_H__
#define __CGLA_VEC2I_H__

#include "ArithVec.h"

namespace CGLA {
	class Vec2f;

	/** \brief 2D Integer vector. */
	
	class Vec2i: public ArithVec<int,Vec2i,2>
	{
	public:
		
		/// Construct 0 vector
		Vec2i() {}

		/// Construct 2D int vector
		Vec2i(int _a,int _b): ArithVec<int,Vec2i,2>(_a,_b) {}

		/// Construct a 2D integer vector with two identical coordinates.
		explicit Vec2i(int a): ArithVec<int,Vec2i,2>(a,a) {}

		/// Convert from 2D float vector
		explicit Vec2i(const Vec2f& v);
  
	};

}
#endif
