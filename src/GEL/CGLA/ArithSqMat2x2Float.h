/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithSqMat2x2Float.h
 * abstract 2x2 floating point matrix class
 */

#ifndef CGLA_ARITHSQMAT2X2FLOAT_H
#define CGLA_ARITHSQMAT2X2FLOAT_H

#include <GEL/CGLA/ExceptionStandard.h>
#include <GEL/CGLA/ArithSqMatFloat.h>


namespace CGLA
{
CGLA_DERIVEEXCEPTION(Mat2x2fException);

/** \brief Two by two float matrix template.

  This class template is useful for various vector transformations
  in the plane. */
template <class V, class M>
class ArithSqMat2x2Float : public ArithSqMatFloat<V, M, 2> {
public:
    /// Vector type
    using VectorType = V;

    /// The type of a matrix element
    using ScalarType = typename V::ScalarType;

public:
    /// Construct a Mat2x2f from two Vec2f vectors.
    constexpr ArithSqMat2x2Float(V a, V b):
        ArithSqMatFloat<V, M, 2>(a, b) {}

    /// Construct a Mat2x2f from four scalars.
    constexpr ArithSqMat2x2Float(ScalarType a, ScalarType b,
                                 ScalarType c, ScalarType d):
        ArithSqMatFloat<V, M, 2>(V(a, b), V(c, d)) {}

    /// Construct the NAN matrix
    constexpr ArithSqMat2x2Float() = default;

    /// Construct a matrix from a single scalar value.
    constexpr explicit ArithSqMat2x2Float(ScalarType a):
        ArithSqMatFloat<V, M, 2>(a) {}
};

/// @name Non-member operations
/// @related ArithSqMat2x2Float
/// @{

/** Compute the determinant of a Mat2x2f. This function is faster than
    the generic determinant function for ArithSqMat */
template <class V, class M>
constexpr typename ArithSqMat2x2Float<V, M>::ScalarType
determinant(const ArithSqMat2x2Float<V, M>& m)
{
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

template <class V, class M>
constexpr M invert(const ArithSqMat2x2Float<V, M>& m)
{
    typename M::ScalarType det = determinant(m);
    if (!is_tiny(fabs(det))) {
        return M(m[1][1] / det, -m[0][1] / det, -m[1][0] / det, m[0][0] / det);
    }
    throw Mat2x2fException("Cannot invert matrix");
}

/// @}

}
#endif
