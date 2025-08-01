/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithVec4Float.h
 * @brief Abstract 4D floating point vector class
 */

#ifndef CGLA_ARITHVEC4FLOAT_H
#define CGLA_ARITHVEC4FLOAT_H

#include <GEL/CGLA/ArithVecFloat.h>

namespace CGLA
{
/** \brief A four-dimensional floating point vector template.

        This class is also used (via typedef) for
        homogeneous vectors.
*/

template <class T, class V>
class ArithVec4Float : public ArithVecFloat<T, V, 4> {
public:
    /// Construct a (0,0,0,0) homogenous Vector
    constexpr ArithVec4Float() = default;

    /// Construct a 4D vector
    constexpr ArithVec4Float(T a, T b, T c, T d):
        ArithVecFloat<T, V, 4>(a, b, c, d) {}

    /// Divide all coordinates by the fourth coordinate
    constexpr void de_homogenize()
    {
        float k = 1.0 / (*this)[3];
        (*this)[0] = (*this)[0] * k;
        (*this)[1] = (*this)[1] * k;
        (*this)[2] = (*this)[2] * k;
        (*this)[3] = 1.0;
    }
};
}
#endif
