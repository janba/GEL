/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <algorithm>

#include "Vec3f.h"
#include "Vec3d.h"
#include "Vec4f.h"

using namespace std;

namespace CGLA 
{

  Vec3f::Vec3f(const Vec3d& v):
    ArithVec3Float<float,Vec3f>(static_cast<float>(v[0]), 
				static_cast<float>(v[1]), 
				static_cast<float>(v[2])) {}
  
  Vec3f::Vec3f(const Vec4f& v):
    ArithVec3Float<float,Vec3f>(v[0],v[1],v[2]) {}

}
