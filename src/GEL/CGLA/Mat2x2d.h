/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pd
 * ----------------------------------------------------------------------- */

/** @file Mat2x2d.h
 * @brief 2x2 double matrix class
 */

#ifndef CGLA_MAT2X2D_H
#define CGLA_MAT2X2D_H

#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/ArithSqMat2x2Float.h>

namespace CGLA
{
/** \brief Two by two double matrix.

  This class is useful for various
  vector transformations in the plane. */
class Mat2x2d : public ArithSqMat2x2Float<Vec2d, Mat2x2d> {
public:
    /// Construct a Mat2x2d from two Vec2d vectors.
    constexpr Mat2x2d(Vec2d a, Vec2d b): ArithSqMat2x2Float(a, b) {}

    /// Construct a Mat2x2f from four scalars.
    constexpr Mat2x2d(double a, double b, double c, double d):
        ArithSqMat2x2Float(Vec2d(a, b), Vec2d(c, d)) {}

    /// Construct the NAN matrix
    constexpr Mat2x2d() = default;

    /// Construct a Mat2x2d from a single scalar
    constexpr explicit Mat2x2d(double a):
        ArithSqMat2x2Float<Vec2d, Mat2x2d>(a) {}
};
}
#endif
