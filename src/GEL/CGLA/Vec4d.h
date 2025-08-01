/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec4d.h
 * @brief 4D double vector class.
 */

#ifndef CGLA_VEC4D_H
#define CGLA_VEC4D_H

#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/ArithVec4Float.h>

namespace CGLA
{
/** \brief A four dimensional floating point vector.

        This class is also used (via typedef) for
        homogeneous vectors.
*/
class Vec4d : public ArithVec4Float<double, Vec4d> {
public:
    /// Construct a (0,0,0,0) homogenous Vector
    constexpr Vec4d(): ArithVec4Float<double, Vec4d>(0.0, 0.0, 0.0, 0.0) {}

    /// Construct a (0,0,0,0) homogenous Vector
    constexpr explicit Vec4d(double _a): ArithVec4Float<double, Vec4d>(_a, _a, _a, _a) {}

    /// Construct a 4D vector
    constexpr Vec4d(double _a, double _b, double _c, double _d):
        ArithVec4Float<double, Vec4d>(_a, _b, _c, _d) {}

    /// Construct a homogenous vector from a non-homogenous.
    constexpr explicit Vec4d(const Vec3d& v, double _d):
        ArithVec4Float<double, Vec4d>(v[0], v[1], v[2], _d) {}
};

class Mat4x4d;

template <>
struct VecT_to_MatT<Vec4d> {
    using MatT = Mat4x4d;
};
}
#endif
