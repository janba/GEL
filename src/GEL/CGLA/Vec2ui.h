/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief 2D unsigned integer vector class.
 * ----------------------------------------------------------------------- */

/** @file Vec2ui.h
 * @brief 2D unsigned integer vector class.
 */

#ifndef CGLA_VEC2UI_H
#define CGLA_VEC2UI_H

#include <GEL/CGLA/ArithVec.h>

namespace CGLA 
{
	class Vec2f;

	/** \brief 2D Integer vector. */
	
	class Vec2ui: public ArithVec<unsigned int,Vec2ui,2>
	{
	public:
		
		/// Construct 0 vector
		constexpr Vec2ui() = default;

		/// Construct 2D int vector
		constexpr Vec2ui(unsigned int _a)
		  : ArithVec<unsigned int,Vec2ui,2>(_a,_a) 
		{}

		/// Construct 2D int vector
		constexpr Vec2ui(unsigned int _a, unsigned int _b)
		  : ArithVec<unsigned int,Vec2ui,2>(_a,_b) 
		{}

		/// Convert from 2D float vector
		constexpr explicit Vec2ui(const Vec2f& v);
  
	};
}
#endif
