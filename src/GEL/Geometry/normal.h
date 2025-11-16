#ifndef GEOMETRY_NORMAL_H
#define GEOMETRY_NORMAL_H
#pragma once

#include <GEL/CGLA/Mat.h>

namespace Geometry
{

CGLA::Vec3d smallestEigenVector(const CGLA::Mat3x3d& matrix);

/// Given an input range of coordinates, estimate a normal vector
template <std::ranges::input_range Range>
CGLA::Vec3d estimateNormal(Range&& neighbors)
{
    CGLA::Vec3d centroid(0.0f);
    for (const auto& point : neighbors) {
        centroid += point;
    }
    centroid /= static_cast<double>(neighbors.size());

    CGLA::Mat3x3d covariance(0.0f);
    for (const auto& point : neighbors) {
        CGLA::Vec3d diff = point - centroid;
        covariance += outer_product(diff, diff);
    }
    covariance /= static_cast<double>(neighbors.size());

    const CGLA::Vec3d normal = smallestEigenVector(covariance);

    return normalize(normal); // Normalize to ensure unit length
}
}

#endif
