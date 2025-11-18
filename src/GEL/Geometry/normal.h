#ifndef GEOMETRY_NORMAL_H
#define GEOMETRY_NORMAL_H
#pragma once

#include <GEL/Util/Assert.h>
#include <GEL/CGLA/Mat.h>
#include <GEL/CGLA/eigensolution.h>

namespace Geometry
{

CGLA::Vec3d smallestEigenVector(const CGLA::Mat3x3d& matrix);

/// Given an input range of coordinates, estimate a normal vector
template <std::ranges::input_range Range>
CGLA::Vec3d estimateNormal(Range&& neighbors, double radius)
{
    CGLA::Vec3d centroid(0.0f);
    for (const auto& point : neighbors) {
        centroid += point;
    }
    centroid /= static_cast<double>(neighbors.size());

    CGLA::Mat3x3d covariance(0.0f);
    for (const auto& point : neighbors) {
        CGLA::Vec3d diff = point - centroid;
        covariance += exp(-2*CGLA::sqr_length(diff)/ CGLA::sqr(radius)) * CGLA::outer_product(diff, diff);
    }
    CGLA::Mat3x3d eigenvectors(0);
    CGLA::Mat3x3d eigenvalues(0);
    // Call GEL's power_eigensolution on the inverse matrix
    auto n = CGLA::power_eigensolution<CGLA::Mat3x3d>(covariance, eigenvectors, eigenvalues, 2);
    GEL_ASSERT(n == 2);
    return CGLA::normalize(CGLA::cross(eigenvectors[0], eigenvectors[1]));
}

}

#endif
