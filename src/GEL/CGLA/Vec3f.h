/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3f.h
 * @brief 3D float vector class.
 */

#ifndef CGLA_VEC3F_H
#define CGLA_VEC3F_H

#include <GEL/CGLA/ArithVec3Float.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3usi.h>

namespace CGLA
{
class Vec3d;
class Vec4f;

/** \brief 3D float vector.

        Class Vec3f is the vector typically used in 3D computer graphics.
        The class has many constructors since we may need to convert from
        other vector types. Most of these are explicit to avoid automatic
        conversion.
*/
class Vec3f : public ArithVec3Float<float, Vec3f> {
public:
    /// Construct 0 vector.
    constexpr Vec3f() = default;

    /// Construct a 3D float vector.
    constexpr Vec3f(float a, float b, float c):
        ArithVec3Float(a, b, c) {}

    /// Construct a vector with 3 identical coordinates.
    constexpr explicit Vec3f(float a):
        ArithVec3Float(a, a, a) {}

    /// Construct from a 3D int vector
    constexpr explicit Vec3f(const Vec3i& v):
        ArithVec3Float(v[0], v[1], v[2]) {}

    /// Construct from a 3D unsigned int vector.
    constexpr explicit Vec3f(const Vec3usi& v):
        ArithVec3Float(v[0], v[1], v[2]) {}

    /// Construct from a 3D double vector.
    explicit Vec3f(const Vec3d&);

    /// Construct from a 4D float vector (skipping the last value)
    explicit Vec3f(const Vec4f&);
};


class Mat3x3f;

template <>
struct VecT_to_MatT<Vec3f> {
    using MatT = Mat3x3f;
};
}
#endif
