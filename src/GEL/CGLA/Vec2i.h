/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec2i.h
 * @brief 2D int vector class.
 */

#ifndef CGLA_VEC2I_H
#define CGLA_VEC2I_H

#include <GEL/CGLA/ArithVec.h>

namespace CGLA
{
class Vec2f;

/// @brief 2D Integer vector.
class Vec2i : public ArithVec<int, Vec2i, 2> {
public:
    /// Construct 0 vector
    constexpr Vec2i() = default;

    /// Construct 2D int vector
    constexpr Vec2i(int _a, int _b): ArithVec(_a, _b) {}

    /// Construct a 2D integer vector with two identical coordinates.
    constexpr explicit Vec2i(int a): ArithVec(a, a) {}

    /// Convert from 2D float vector
    explicit Vec2i(const Vec2f& v);
};
}
#endif
