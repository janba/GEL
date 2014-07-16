/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3uc.h
 * @brief 3D unsigned char vector.
 */

#ifndef __CGLA_VEC3UC_H__
#define __CGLA_VEC3UC_H__

#include "Vec3i.h"

namespace CGLA {
	typedef unsigned char UChar;

	/** \brief 3D unsigned char vector. */
	class Vec3uc: public ArithVec3Int<UChar,Vec3uc>
	{

	public:
		
		/// Construct 0 vector
		Vec3uc() {}

		/// Construct 3D uchar vector
		Vec3uc(UChar a, UChar b, UChar c): 
			ArithVec3Int<UChar,Vec3uc>(a,b,c) {}

		/// Convert from int vector. 
		explicit Vec3uc(const Vec3i& v): 
			ArithVec3Int<UChar,Vec3uc>(v[0]&0xff, v[1]&0xff, v[2]&0xff) {}
	};


}
#endif

