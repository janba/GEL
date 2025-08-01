/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3d.h
 * @brief 3D double vector class.
 */

#ifndef CGLA_VEC3D_H
#define CGLA_VEC3D_H

#include <GEL/CGLA/ArithVec3Float.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3usi.h>
#include <GEL/CGLA/Vec3f.h>


namespace CGLA
{
/** \brief A 3D double vector.
Useful for high precision arithmetic. */
class Vec3d : public ArithVec3Float<double, Vec3d> {
public:
    /// Construct 0 vector
    constexpr Vec3d() = default;

    /// Construct vector
    constexpr Vec3d(double a, double b, double c): ArithVec3Float<double, Vec3d>(a, b, c) {}

    /// Construct vector where all coords = a
    constexpr explicit Vec3d(double a):
        ArithVec3Float<double, Vec3d>(a, a, a) {}

    /// Convert from int vector
    constexpr explicit Vec3d(const Vec3i& v):
        ArithVec3Float<double, Vec3d>(v[0], v[1], v[2]) {}

    /// Construct from a 3D unsigned int vector.
    constexpr explicit Vec3d(const Vec3usi& v):
        ArithVec3Float<double, Vec3d>(v[0], v[1], v[2]) {}

    /// Convert from float vector
    constexpr explicit Vec3d(const Vec3f& v):
        ArithVec3Float<double, Vec3d>(v[0], v[1], v[2]) {}
};

class Mat3x3d;

template <>
struct VecT_to_MatT<Vec3d> {
    using MatT = Mat3x3d;
};
}
#endif
