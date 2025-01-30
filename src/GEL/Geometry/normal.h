#ifndef Normal_h
#define Normal_h
#pragma once

#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/eigensolution.h>
#include <GEL/CGLA/Mat3x3d.h>

using namespace CGLA;

namespace Geometry {
	Vec3d estimateNormal(const std::vector<Vec3d>& neighbors);

	Vec3f smallestEigenVector(const Mat3x3f& matrix);
}

#endif