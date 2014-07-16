/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file   Vec4uc.h
 * @brief 4D unsigned char vector
 */

#ifndef __CGLA_VEC4UC_H__
#define __CGLA_VEC4UC_H__

#include "Vec4f.h"

namespace CGLA {
	typedef unsigned char UChar;

	/** \brief 4D unsigned char vector. */
	class Vec4uc: public ArithVec<UChar,Vec4uc,4>
	{

	public:
		
		/// Construct 0 vector
		Vec4uc() {}

		/// Construct 0 vector
		Vec4uc(unsigned char a): ArithVec<UChar,Vec4uc,4>(a,a,a,a) {}

		/// Construct 4D uchar vector
		Vec4uc(UChar _a, UChar _b, UChar _c,UChar _d): 
			ArithVec<UChar,Vec4uc,4>(_a,_b,_c,_d) {}

		/// Convert from float vector. 
		explicit Vec4uc(const Vec4f& v): 
			ArithVec<UChar,Vec4uc,4>(UChar(v[0]), UChar(v[1]), 
															 UChar(v[2]), UChar(v[3])) {}

		operator Vec4f() const
		{
			return  Vec4f((*this)[0],(*this)[1],(*this)[2],(*this)[3]);
		}

	};


}
#endif

