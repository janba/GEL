#ifndef GEOMETRY_NORMAL_H
#define GEOMETRY_NORMAL_H
#pragma once

#include <GEL/CGLA/Mat.h>
#include <GEL/CGLA/eigensolution.h>
#include <GEL/CGLA/Vec.h>
#include <GEL/CGLA/ls_solve.h>

namespace Geometry
{
#include <array>
#include <cmath>
#include <iostream>

/// Given an input range of coordinates, estimate a normal vector
template <std::ranges::input_range Range>
CGLA::Vec3d estimateNormal(Range&& neighbors, double radius)
{
    CGLA::Vec3d centroid(0.0f);
    for (const auto& point : neighbors)
        centroid += point;
    centroid /= static_cast<double>(neighbors.size());

    CGLA::Mat3x3d covariance(0.0f);
    for (const auto& point : neighbors) {
        CGLA::Vec3d diff = point - centroid;
        covariance += CGLA::outer_product(diff, diff);// * exp(-4*CGLA::sqr_length(diff)/(radius*radius));
    }

    CGLA::Mat3x3d eigenvectors(0), eigenvalues(0);
    CGLA::power_eigensolution<CGLA::Mat3x3d>(covariance, eigenvectors, eigenvalues, 2);
    CGLA::Vec3d norm = CGLA::normalize(CGLA::cross(eigenvectors[0], eigenvectors[1]));
    return norm;
    // CGLA::Vec3d T,B;
    // CGLA::orthogonal(norm, T, B);

    // CGLA::Mat3x3d M(0.0);
    // CGLA::Vec3d b(0.0);
    // for (auto& point : neighbors) {
    //     CGLA::Vec3d diff = point - centroid;
    //     double u = CGLA::dot(diff, T) / radius;
    //     double v = CGLA::dot(diff, B) / radius;
    //     double z = CGLA::dot(diff, norm) / radius;
    //     M += CGLA::Mat3x3d(
    //         CGLA::Vec3d(1, u, v),
    //         CGLA::Vec3d(u, u * u, u * v),
    //         CGLA::Vec3d(v, v * u, v * v));
    //     b[0] += z;
    //     b[1] += u * z;
    //     b[2] += v * z;
    // }
    // CGLA::Vec3d x = CGLA::lin_solve(M, b);
    // return CGLA::normalize(-x[1] * T - x[2] * B + norm);

}
}

#endif
