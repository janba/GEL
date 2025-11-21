/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithVecFloat.h
 * @brief Abstract 2D floating point vector class
 */

#ifndef CGLA_ARITHVECFLOAT_H
#define CGLA_ARITHVECFLOAT_H

#include <GEL/CGLA/ArithVec.h>

namespace CGLA
{
template <class T, class V, unsigned int N>
class ArithVecFloat : public ArithVec<T, V, N> {
public:
    constexpr ArithVecFloat()
    {
        std::fill_n(this->data.begin(), N, CGLA_INIT_VALUE);
    }

    constexpr explicit ArithVecFloat(T a):
        ArithVec<T, V, N>(a) {}

    constexpr ArithVecFloat(T a, T b):
        ArithVec<T, V, N>(a, b) {}

    constexpr ArithVecFloat(T a, T b, T c):
        ArithVec<T, V, N>(a, b, c) {}

    constexpr ArithVecFloat(T a, T b, T c, T d):
        ArithVec<T, V, N>(a, b, c, d) {}

    /// Compute Euclidean length.
    [[nodiscard]]
    constexpr T length() const
    {
        return std::sqrt(sqr_length(*this));
    }

    /// Normalize vector.
    constexpr void normalize()
    {
        (*this) /= this->length();
    }

    /// Conditionally normalize vector. The condition being that the vector has non-zero length
    constexpr void cond_normalize()
    {
        T sql = sqr_length(*this);
        if (sql > 0)
            (*this) /= std::sqrt(sql);
    }
};

/// @name Non-member operations
/// @related ArithVecFloat
/// @{

/// Returns length of vector
template <class T, class V, unsigned int N>
[[nodiscard]]
constexpr T length(const ArithVecFloat<T, V, N>& v)
{
    return v.length();
}


/// Returns normalized vector
template <class T, class V, unsigned int N>
[[nodiscard]]
constexpr V normalize(const ArithVecFloat<T, V, N>& v)
{
    return v / v.length();
}

/// Returns normalized vector if the vector has non-zero length - otherwise the 0 vector.
template <class T, class V, unsigned int N>
[[nodiscard]]
constexpr V cond_normalize(const ArithVecFloat<T, V, N>& v)
{
    T sql = sqr_length(v);
    if (sql > 0)
        return v / std::sqrt(sql);
    return v * 1.0;
}

/// Returns the linear interpolation of v1 and v2 given t in [0, 1].
/// Values of t out of this range will still return a valid value but probably does not do what you want
template <class T, class V, unsigned int N>
constexpr V lerp(const ArithVecFloat<T, V, N>& v1, const ArithVecFloat<T, V, N>& v2, double t)
{
    return v1 * (1.0 - t) + v2 * t;
}

/// @}

/** The template below is used to map vector types to matrix types. In each of the floating point vector classes
 the template is specialized such that when the template argument is that particular vector class, the appropriate
 matrix class is represented by MatT. For instance, VecT_to_MatT<Vec3d>::MatT is the type Mat3x3d */
/// @related ArithVecFloat
template <typename V>
struct VecT_to_MatT {
    using MatT = void;
};
}

#endif
