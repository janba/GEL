/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithSqMat3x3Float.h
 * @brief Abstract 3x3 floating point matrix class
 */

#ifndef CGLA_ARITHSQMAT3X3FLOAT_H
#define CGLA_ARITHSQMAT3X3FLOAT_H

#include <GEL/CGLA/ExceptionStandard.h>
#include <GEL/CGLA/CGLA-util.h>
#include <GEL/CGLA/ArithSqMatFloat.h>

namespace CGLA
{
CGLA_DERIVEEXCEPTION(Mat3x3fException)

CGLA_DERIVEEXCEPTION(Mat3x3fSingular)

/** \brief 3 by 3 float matrix template.

    This class template will typically be used for rotation or
    scaling matrices for 3D vectors. */
template <class V, class M>
class ArithSqMat3x3Float : public ArithSqMatFloat<V, M, 3> {
public:
    /// Vector type
    using VectorType = V;

    /// The type of matrix element
    using ScalarType = typename V::ScalarType;

public:
    /// Construct matrix from 3 Vec3f vectors.
    constexpr ArithSqMat3x3Float(V _a, V _b, V _c):
        ArithSqMatFloat<V, M, 3>(_a, _b, _c) {}

    /// Construct the 0 matrix
    constexpr ArithSqMat3x3Float() = default;

    /// Construct a matrix from a single scalar value.
    constexpr explicit ArithSqMat3x3Float(ScalarType a):
        ArithSqMatFloat<V, M, 3>(a) {}
};

/// @name Non-member operations
/// @related ArithSqMat3x3Float
/// @{

/// Invert 3x3 matrix
template <class V, class M>
constexpr M invert(const ArithSqMat3x3Float<V, M>& _a)
{
    // From Graphics Gems IV, Jean-Francois Doue, C++ Vector and Matrix Algebra Routines.
    // From the EULA "Using the code is permitted in any program, product, or library, non-commercial
    // or commercial. Giving credit is not required, though is a nice gesture"
    M a(_a[0], _a[1], _a[2]);
    M b;
    b.identity();

    for (unsigned int j = 0; j < 3; j++) {
        unsigned int i1 = j; // Row with largest pivot candidate
        for (unsigned int i = j + 1; i < 3; i++)
            if (fabs(a[i][j]) > fabs(a[i1][j]))
                i1 = i;

        // Swap rows i1 and j in a and b to put pivot on diagonal
        V a_tmp = a[i1];
        a[i1] = a[j];
        a[j] = a_tmp;

        V b_tmp = b[i1];
        b[i1] = b[j];
        b[j] = b_tmp;

        // Scale row j to have a unit diagonal
        if (a[j][j] == 0.0)
            throw Mat3x3fSingular("Tried to invert Singular matrix");

        b[j] /= a[j][j];
        a[j] /= a[j][j];

        // Eliminate off-diagonal elems in col j of a, doing identical ops to b
        for (unsigned int i = 0; i < 3; i++)
            if (i != j) {
                b[i] -= a[i][j] * b[j];
                a[i] -= a[i][j] * a[j];
            }
    }
    return b;
}

/** Compute determinant. There is a more generic function for
    computing determinants of square matrices (ArithSqMat). This one
    is faster but works only on Mat3x3f */
template <class V, class M>
constexpr typename ArithSqMat3x3Float<V, M>::ScalarType
determinant(const ArithSqMat3x3Float<V, M>& m)
{
    return
        m[0][0] * m[1][1] * m[2][2] +
        m[0][1] * m[1][2] * m[2][0] +
        m[0][2] * m[1][0] * m[2][1] -
        m[0][2] * m[1][1] * m[2][0] -
        m[0][0] * m[1][2] * m[2][1] -
        m[0][1] * m[1][0] * m[2][2];
}

/// @}

}
#endif
