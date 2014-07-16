/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3Hf.h
 * @brief Homogeneous 3D vector class - really a 4D vector.
 */

#ifndef __CGLA_VEC3HF_H__
#define __CGLA_VEC3HF_H__

#include "Vec4f.h"

namespace CGLA {

	/** A 3D homogeneous vector is simply a four D vector.
			I find this simpler than a special class for homogeneous
			vectors. */
	typedef Vec4f Vec3Hf;

}
#endif

