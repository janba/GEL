#ifndef GEOMETRY_NORMAL_H
#define GEOMETRY_NORMAL_H
#pragma once

#include <GEL/Util/Assert.h>
#include <GEL/CGLA/Mat.h>
#include <GEL/CGLA/eigensolution.h>

namespace Geometry
{
#include <array>
#include <cmath>
#include <iostream>
#include <utility>

using EigenResult = std::pair<CGLA::Mat3x3d, CGLA::Mat3x3d>;

inline EigenResult jacobiEigenDecomposition(const CGLA::Mat3x3d & A, int maxIter = 50, double tol = 1e-12) {
    // Initialize eigenvectors as identity
    CGLA::Mat3x3d V = CGLA::identity_Mat3x3d();
    CGLA::Mat3x3d L = A;

    for (int iter = 0; iter < maxIter; ++iter) {
        // Find largest off-diagonal element
        int p = 0, q = 1;
        double maxVal = std::fabs(L[0][1]);
        for (int i = 0; i < 3; ++i) {
            for (int j = i+1; j < 3; ++j) {
                if (std::fabs(L[i][j]) > maxVal) {
                    maxVal = std::fabs(L[i][j]);
                    p = i; q = j;
                }
            }
        }
        if (maxVal < tol) break; // Converged

        // Compute rotation
        double theta = 0.5 * std::atan2(2*L[p][q], L[q][q] - L[p][p]);
        double c = std::cos(theta);
        double s = std::sin(theta);

        // Rotate A
        for (int k = 0; k < 3; ++k) {
            double Apk = L[p][k], Aqk = L[q][k];
            L[p][k] = c*Apk - s*Aqk;
            L[q][k] = s*Apk + c*Aqk;
        }
        for (int k = 0; k < 3; ++k) {
            double Akp = L[k][p], Akq = L[k][q];
            L[k][p] = c*Akp - s*Akq;
            L[k][q] = s*Akp + c*Akq;
        }

        // Rotate V
        for (int k = 0; k < 3; ++k) {
            double Vkp = V[k][p], Vkq = V[k][q];
            V[k][p] = c*Vkp - s*Vkq;
            V[k][q] = s*Vkp + c*Vkq;
        }
    }
    V = CGLA::transpose(V);
    return std::make_pair(V, L);
}




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
        covariance += CGLA::outer_product(diff, diff);
    }
    // auto [Q, L] = jacobiEigenDecomposition(covariance);
    // return CGLA::normalize(Q[0]);
    CGLA::Mat3x3d eigenvectors(0);
    CGLA::Mat3x3d eigenvalues(0);
    // Call GEL's power_eigensolution on the inverse matrix
    auto n = CGLA::power_eigensolution<CGLA::Mat3x3d>(covariance, eigenvectors, eigenvalues, 2);
    GEL_ASSERT(n == 2);
    return CGLA::normalize(CGLA::cross(eigenvectors[0], eigenvectors[1]));

    // CGLA::Vec3d n(1.0, 1.0, 1.0);
    // CGLA::Vec3d n_old = n;
    // for (int i = 0; i<1000; ++i) {
    //     for (int c0_idx = 0; c0_idx < 3; ++c0_idx) {
    //         int c1_idx = (i + 1) % 3;
    //         int c2_idx = (i + 2) % 3;
    //         for (int p_idx = 0; p_idx < neighbors.size(); ++p_idx) {
    //             CGLA::Vec3d p = neighbors[p_idx] - centroid;
    //             n[c0_idx] -= (n[c1_idx] * p[c1_idx] + n[c2_idx] * p[c2_idx]) / p[c0_idx];
    //         }
    //         n[c0_idx] /= static_cast<double>(neighbors.size());
    //     }
    //     n = CGLA::normalize(n);
    //     if (dot(n, n_old) > 0.9999)
    //         break;
    //     n_old = n;
    // }
    // return n;

}

// n.x * p_i.x + n.y * p_i.y + n.z * p_i.z = 0
// n.x = - (n.y * p_i.y + n.z * p_i.z) / p_i.x
}

#endif
