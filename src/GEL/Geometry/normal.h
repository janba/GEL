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

inline bool solve6x6(double A[6][6], double B[6], double X[6]) {
    // Augmented matrix
    float M[6][7];
    for (int r = 0; r < 6; r++) {
        for (int c = 0; c < 6; c++) M[r][c] = A[r][c];
        M[r][6] = B[r];
    }

    // Gaussian elimination with partial pivoting
    for (int k = 0; k < 6; k++)
    {
        // Pivot
        int pivot = k;
        float maxv = fabs(M[k][k]);
        for (int r = k + 1; r < 6; r++)
            if (fabs(M[r][k]) > maxv) { maxv = fabs(M[r][k]); pivot = r; }

        if (maxv < 1e-12f) return false;

        if (pivot != k)
            for (int c = 0; c < 7; c++) std::swap(M[k][c], M[pivot][c]);

        float div = M[k][k];
        for (int c = k; c < 7; c++) M[k][c] /= div;

        for (int r = 0; r < 6; r++)
        {
            if (r == k) continue;
            float f = M[r][k];
            for (int c = k; c < 7; c++)
                M[r][c] -= f * M[k][c];
        }
    }

    for (int i = 0; i < 6; i++)
        X[i] = M[i][6];

    return true;
}

using namespace CGLA;

template <std::ranges::input_range Range>
Vec3d estimateNormalJet(Range&& neighbors, double radius) {
    CGLA::Vec3d centroid(0.0f);
    for (const auto& point : neighbors)
        centroid += point;
    centroid /= static_cast<double>(neighbors.size());

    Vec3d nPCA = estimateNormal(neighbors, radius); 

    // Build local frame
    Vec3d tmp = fabs(nPCA[0]) < 0.9f ? Vec3d(1, 0, 0) : Vec3d(0, 1, 0);
    Vec3d u = cross(nPCA, tmp); 
    u.normalize();
    Vec3d v = cross(nPCA, u); 
    v.normalize();

    // Prepare arrays
    double A[6][6] = { 0 };
    double B[6] = { 0 };

    for (const auto& point : neighbors)
    {
        Vec3d d = point - centroid;

        double uu = dot(d, u);
        double vv = dot(d, v);
        double ww = dot(d, nPCA);

        double row[6] = { 1, uu, vv, uu * uu, uu * vv, vv * vv };

        for (int r = 0; r < 6; r++)
        {
            for (int c = 0; c < 6; c++) A[r][c] += row[r] * row[c];
            B[r] += row[r] * ww;
        }
    }

    double X[6];
    if (!solve6x6(A, B, X))
        return nPCA;

    // Get local normal
    double a1 = X[1], a2 = X[2];
    Vec3d n_local(-a1, -a2, 1.);
    if (sqr_length(n_local) < 1e-10f)
        return nPCA;
    n_local.normalize();

    Vec3d output = Vec3d(n_local[0] * u +
        n_local[1] * v +
        n_local[2] * nPCA);
    output.normalize();
    // Back to world coords
    return output;
}

}

#endif
