/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat2x3f.h
 * @brief 2x3 float matrix class
 */

#ifndef CGLA_MAT2X3F_H
#define CGLA_MAT2X3F_H

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/ArithMatFloat.h>

namespace CGLA
{
// /**  \brief 2x3 float matrix class.
//
//      This class is useful for projecting a vector from 3D space to 2D.
// */
// class Mat2x3f : public ArithMatFloat<Vec2f, Vec3f, Mat2x3f, 2> {
// public:
//     /// Construct Mat2x3f from two Vec3f vectors (vectors become rows)
//     constexpr Mat2x3f(const Vec3f& _a, const Vec3f& _b):
//         ArithMatFloat<Vec2f, Vec3f, Mat2x3f, 2>(_a, _b) {}
//
//     /// Construct NAN matrix.
//     constexpr Mat2x3f() = default;
// };
//
//
// /**  \brief 3x2 float matrix class.
//
//      This class is useful for going from plane to 3D coordinates.
// */
// class Mat3x2f : public ArithMatFloat<Vec3f, Vec2f, Mat3x2f, 3> {
// public:
//     /** Construct matrix from three Vec2f vectors which become the
//     rows of the matrix. */
//     constexpr Mat3x2f(const Vec2f& _a, const Vec2f& _b, const Vec2f& _c):
//         ArithMatFloat<Vec3f, Vec2f, Mat3x2f, 3>(_a, _b, _c) {}
//
//     /// Construct NAN matrix.
//     constexpr Mat3x2f() {}
// };
}
#endif
