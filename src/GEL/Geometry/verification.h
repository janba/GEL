/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file verification.h
 * @brief Test implementation of distance to triangle. Only for reference.
 */
#ifndef __GEOMETRY_VERIFICATION_H
#define __GEOMETRY_VERIFICATION_H

#include "Triangle.h"

namespace Geometry
{
  float SqrDistance (const CGLA::Vec3f& rkPoint,const Triangle& rkTri);
}
#endif
