/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract 2D floating point vector class
 * ----------------------------------------------------------------------- */

/** @file ArithVec2Float.h
 * @brief Abstract 2D floating point vector class
 */

#ifndef CGLA__ARITHVEC2FLOAT_H
#define CGLA__ARITHVEC2FLOAT_H

#include <GEL/CGLA/ArithVecFloat.h>

namespace CGLA
{
template <class T, class V>
class ArithVec2Float : public ArithVecFloat<T, V, 2> {
public:
    /// Construct a 2D float vector.
    constexpr ArithVec2Float(T a, T b): ArithVecFloat<T, V, 2>(a, b) {}

    /// Construct a 2D float vector.
    constexpr ArithVec2Float() = default;
};

/// @name Non-member operations
/// @related ArithVec2Float
/// @{

/// Rotates vector 90 degrees to obtain orthogonal vector
template <class T, class V>
constexpr V orthogonal(const ArithVec2Float<T, V>& v)
{
    return V(-v[1], v[0]);
}

// Computes (scalar) cross product from two vectors
template <class T, class V>
constexpr T cross(const ArithVec2Float<T, V>& a,
                  const ArithVec2Float<T, V>& b)
{
    return a[0] * b[1] - a[1] * b[0];
}

/// The two last (scalar) arguments are the linear combination of
/// the two first arguments (vectors) which produces the third argument.
template<class T, class V>
bool linear_combine(const ArithVec2Float<T,V>& a,
                                        const ArithVec2Float<T,V>& b,
                                        const ArithVec2Float<T,V>& c,
                                        T&,
                                        T&);

/// @}
}

#endif
