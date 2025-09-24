//
// Created by Cem Akarsubasi on 8/4/25.
//

#ifndef GEL_MAT_H
#define GEL_MAT_H

#include <GEL/CGLA/Definitions.h>
#include <GEL/CGLA/ArithSqMat2x2Float.h>
#include <GEL/CGLA/ArithSqMat3x3Float.h>
#include <GEL/CGLA/ArithSqMat4x4Float.h>

namespace CGLA
{
/// @brief Two by two double matrix.
///
///  This class is useful for various
///  vector transformations in the plane.
template <std::floating_point Float>
class Mat2x2 final : public ArithSqMat2x2Float<Vec2<Float>, Mat2x2<Float>> {
public:
    /// Construct a Mat2x2 from two Vec2 vectors.
    constexpr Mat2x2(const Vec2<Float>& a, const Vec2<Float>& b) : ArithSqMat2x2Float<Vec2<Float>, Mat2x2d>(a, b) {}

    /// Construct a Mat2x2 from four scalars.
    constexpr Mat2x2(Float a, Float b, Float c, Float d) :
        ArithSqMat2x2Float<Vec2<Float>, Mat2x2>(Vec2(a, b), Vec2(c, d)) {}

    /// Construct the NAN matrix
    constexpr Mat2x2() = default;

    /// Construct a Mat2x2d from a single scalar
    constexpr explicit Mat2x2(Float a) :
        ArithSqMat2x2Float<Vec2<Float>, Mat2x2>(a) {}
};


/// @brief 3 by 3 double matrix.
///
/// This class will typically be used for rotation or
/// scaling matrices for 3D vectors.
template <std::floating_point Float>
class Mat3x3 final : public ArithSqMat3x3Float<Vec3<Float>, Mat3x3<Float>> {
public:
    /// Construct matrix from 3 Vec3d vectors.
    constexpr Mat3x3(const Vec3<Float>& _a, const Vec3<Float>& _b, const Vec3<Float>& _c) :
        ArithSqMat3x3Float<Vec3<Float>, Mat3x3>(_a, _b, _c) {}

    /// Construct the 0 matrix
    constexpr Mat3x3() = default;

    /// Construct a matrix from a single scalar value.
    constexpr explicit Mat3x3(Float a) : ArithSqMat3x3Float<Vec3<Float>, Mat3x3>(a) {}
};


/** \brief 4x4 double matrix.

    This class is useful for transformations such as perspective projections
    or translation where 3x3 matrices do not suffice. */
template <std::floating_point Float>
class Mat4x4 final : public ArithSqMat4x4Float<Vec4<Float>, Mat4x4<Float>> {
public:
    /// Construct a Mat4x4d from four Vec4d vectors
    constexpr Mat4x4(const Vec4<Float>& _a, const Vec4<Float>& _b, const Vec4<Float>& _c, const Vec4<Float>& _d) :
        ArithSqMat4x4Float<Vec4<Float>, Mat4x4>(_a, _b, _c, _d) {}

    /// Construct the nan matrix
    constexpr Mat4x4() = default;

    /// Construct a matrix with identical elements.
    constexpr explicit Mat4x4(Float a) : ArithSqMat4x4Float<Vec4<Float>, Mat4x4>(a) {}
};


///  @brief 2x3 matrix class.
///
///  This class is useful for projecting a vector from 3D space to 2D.
template <std::floating_point Float>
class Mat2x3 final : public ArithMatFloat<Vec2<Float>, Vec3<Float>, Mat2x3<Float>, 2> {
public:
    /// Construct Mat2x3d from two Vec3f vectors (vectors become rows)
    constexpr Mat2x3(const Vec3<Float>& _a, const Vec3<Float>& _b) :
        ArithMatFloat<Vec2<Float>, Vec3<Float>, Mat2x3, 2>(_a, _b) {}

    /// Construct 0 matrix.
    constexpr Mat2x3() = default;
};

/// @brief 3x2 double matrix class.
///
/// This class is useful for going from plane to 3D coordinates.
template <std::floating_point Float>
class Mat3x2 final : public ArithMatFloat<Vec3<Float>, Vec2<Float>, Mat3x2<Float>, 3> {
public:
    /** Construct matrix from three Vec2d vectors which become the
            rows of the matrix. */
    constexpr Mat3x2(const Vec2<Float>& _a, const Vec2<Float>& _b, const Vec2<Float>& _c) :
        ArithMatFloat<Vec3<Float>, Vec2<Float>, Mat3x2, 3>(_a, _b, _c) {}

    /// Construct 0 matrix.
    constexpr Mat3x2() = default;
};
}

#endif //GEL_MAT_H
