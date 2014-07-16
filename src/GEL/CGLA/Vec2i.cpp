/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Vec2i.h"
#include "Vec2f.h"

namespace CGLA {

	Vec2i::Vec2i(const Vec2f& v):
		ArithVec<int,Vec2i,2>(int(floor(v[0])), int(floor(v[1]))) {}

}
