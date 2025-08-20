/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file QEM.h
 * @brief Quadric Error Metrics. To be used when simplifying.
 */

#ifndef GEOMETRY_QEM_H
#define GEOMETRY_QEM_H

#include <GEL/CGLA/Mat.h>

namespace Geometry
{
namespace detail
{
    inline CGLA::Mat3x3d direct_product(const CGLA::Vec3d& v0, const CGLA::Vec3d& v1)
    {
        CGLA::Mat3x3d m;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i][j] = v0[i] * v1[j];
        return m;
    }
}

/// Quadratic Error Metric
class QEM {
    CGLA::Mat3x3d A;
    CGLA::Vec3d b;
    double c;

public:
    QEM() : A(0), b(0), c(0) {}

    /// Construct a QEM from the given position and vector
    QEM(const CGLA::Vec3d& p0, const CGLA::Vec3d& n0, double w = 1.0f) :
        A(detail::direct_product(n0, n0) * w),
        b(-2 * n0 * dot(n0, p0) * w),
        c(dot(p0, n0) * dot(p0, n0) * w) {}


    void operator+=(const QEM& rhs)
    {
        A += rhs.A;
        b += rhs.b;
        c += rhs.c;
    }

    /// Add the metrics together
    QEM operator+(const QEM& rhs) const
    {
        QEM r = *this;
        r += rhs;
        return r;
    }

    /// Calculate the error for a given position on this QEM
    [[nodiscard]]
    double error(const CGLA::Vec3d& p) const
    {
        return CGLA::dot(p, A * p) + CGLA::dot(b, p) + c;
    }

    [[nodiscard]]
    double determinant() const
    {
        return CGLA::determinant(A);
    }

    [[nodiscard]]
    CGLA::Vec3d grad(const CGLA::Vec3d& p) const
    {
        return CGLA::Vec3d(2 * A * p + b);
    }

    /// Calculate the optimal position for a given singularity threshold and vertex coordinate
    [[nodiscard]]
    CGLA::Vec3d opt_pos(double QEM_thresh = 0.5, const CGLA::Vec3d& p0 = CGLA::Vec3d(0.0)) const;
};
}

// FIXME: remove/move this namespace alias
namespace GEO = Geometry;

#endif
