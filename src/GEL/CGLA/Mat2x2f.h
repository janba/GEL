/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat2x2f.h
 * @brief 2x2 float matrix class
 */

#ifndef CGLA_MAT2X2F_H
#define CGLA_MAT2X2F_H

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/ArithSqMat2x2Float.h>


namespace CGLA
{
/** \brief Two by two float matrix.

  This class is useful for various
  vector transformations in the plane. */
class Mat2x2f : public ArithSqMat2x2Float<Vec2f, Mat2x2f> {
public:
    /// Construct a Mat2x2f from two Vec2f vectors.
    constexpr Mat2x2f(Vec2f _a, Vec2f _b): ArithSqMat2x2Float(_a, _b) {}

    /// Construct a Mat2x2f from four scalars.
    constexpr Mat2x2f(float _a, float _b, float _c, float _d):
        ArithSqMat2x2Float(Vec2f(_a, _b), Vec2f(_c, _d)) {}

    /// Construct the NAN matrix
    constexpr Mat2x2f() = default;

    /// Construct a Mat2x2f from a single scalar
    constexpr explicit Mat2x2f(float a): ArithSqMat2x2Float(a) {}
};
}
#endif
