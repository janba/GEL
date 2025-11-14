#ifndef CGLA_Normal_h
#define CGLA_Normal_h
#pragma once

#include <GEL/CGLA/Definitions.h>
#include <vector>

using namespace CGLA;

namespace Geometry {
	Vec3d estimateNormal(const std::vector<Vec3d>& neighbors);

	Vec3f smallestEigenVector(const Mat3x3f& matrix);
}

#endif