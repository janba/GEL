/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec2f.h
 * @brief 2D float vector class.
 */

#ifndef CGLA_VEC2F_H
#define CGLA_VEC2F_H

#include <GEL/CGLA/ArithVec2Float.h>

namespace CGLA
{

/// @brief 2D floating point vector
class Vec2f : public ArithVec2Float<float, Vec2f> {
public:
    constexpr Vec2f() = default;
    constexpr Vec2f(float _a, float _b): ArithVec2Float(_a, _b) {}

    template <class T, class V, unsigned int N>
    constexpr explicit Vec2f(const ArithVec<T, V, N>& v):
        ArithVec2Float<float, Vec2f>(static_cast<float>(v[0]),
                                     static_cast<float>(v[1])) {}

    constexpr explicit Vec2f(float a): ArithVec2Float<float, Vec2f>(a, a) {}
};

class Mat2x2f;

template <>
struct VecT_to_MatT<Vec2f> {
    using MatT = Mat2x2f;
};
}
#endif
