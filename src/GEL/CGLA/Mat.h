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
class Mat2x2 : public ArithSqMat2x2Float<Vec2<Float>, Mat2x2<Float>> {
public:
    /// Construct a Mat2x2d from two Vec2d vectors.
    constexpr Mat2x2(const Vec2<Float>& a, const Vec2<Float>& b): ArithSqMat2x2Float<Vec2<Float>, Mat2x2d>(a, b) {}

    /// Construct a Mat2x2f from four scalars.
    constexpr Mat2x2(Float a, Float b, Float c, Float d):
        ArithSqMat2x2Float<Vec2<Float>, Mat2x2>(Vec2(a, b), Vec2(c, d)) {}

    /// Construct the NAN matrix
    constexpr Mat2x2() = default;

    /// Construct a Mat2x2d from a single scalar
    constexpr explicit Mat2x2(Float a):
        ArithSqMat2x2Float<Vec2<Float>, Mat2x2>(a) {}
};


template <std::floating_point Float>

/// 3 by 3 double matrix.
///
/// This class will typically be used for rotation or
/// scaling matrices for 3D vectors.
class Mat3x3 : public ArithSqMat3x3Float<Vec3<Float>, Mat3x3<Float>> {
public:
    /// Construct matrix from 3 Vec3d vectors.
    constexpr Mat3x3(const Vec3<Float>& _a, const Vec3<Float>& _b, const Vec3<Float>& _c):
        ArithSqMat3x3Float<Vec3<Float>, Mat3x3>(_a, _b, _c) {}

    /// Construct the 0 matrix
    constexpr Mat3x3() = default;

    /// Construct a matrix from a single scalar value.
    constexpr explicit Mat3x3(Float a): ArithSqMat3x3Float<Vec3<Float>, Mat3x3>(a) {}
};


/** \brief 4x4 double matrix.

    This class is useful for transformations such as perspective projections
    or translation where 3x3 matrices do not suffice. */
template <std::floating_point Float>
class Mat4x4 : public ArithSqMat4x4Float<Vec4<Float>, Mat4x4<Float>> {
public:
    /// Construct a Mat4x4d from four Vec4d vectors
    constexpr Mat4x4(const Vec4<Float>& _a, const Vec4<Float>& _b, const Vec4<Float>& _c, const Vec4<Float>& _d):
        ArithSqMat4x4Float<Vec4<Float>, Mat4x4>(_a, _b, _c, _d) {}

    /// Construct the nan matrix
    constexpr Mat4x4() = default;

    /// Construct a matrix with identical elements.
    constexpr explicit Mat4x4(Float a): ArithSqMat4x4Float<Vec4<Float>, Mat4x4>(a) {}
};


///  @brief 2x3 matrix class.
///
///  This class is useful for projecting a vector from 3D space to 2D.
template <std::floating_point Float>
class Mat2x3: public ArithMatFloat<Vec2<Float>, Vec3<Float>, Mat2x3<Float>, 2>
{

public:
    /// Construct Mat2x3d from two Vec3f vectors (vectors become rows)
    constexpr Mat2x3(const Vec3<Float>& _a, const Vec3<Float>& _b):
        ArithMatFloat<Vec2<Float>, Vec3<Float>, Mat2x3, 2> (_a,_b) {}

    /// Construct 0 matrix.
    constexpr Mat2x3() = default;
};

/// @brief 3x2 double matrix class.
///
/// This class is useful for going from plane to 3D coordinates.
template <std::floating_point Float>
class Mat3x2: public ArithMatFloat<Vec3<Float>, Vec2<Float>, Mat3x2<Float>, 3> {
public:

    /** Construct matrix from three Vec2d vectors which become the
            rows of the matrix. */
    constexpr Mat3x2(const Vec2<Float>& _a, const Vec2<Float>& _b, const Vec2<Float>& _c):
        ArithMatFloat<Vec3<Float>, Vec2<Float>, Mat3x2, 3> (_a,_b,_c) {}

    /// Construct 0 matrix.
    constexpr Mat3x2() = default;
};

template <typename Matrix>
constexpr Matrix identity()
{
    Matrix m(0.0);
    static_assert(m.get_v_dim() == m.get_h_dim(), "Matrix must be square");
    for (int i = 0; i < m.get_v_dim(); ++i) {
        m[i][i] = 1.0;
    }
    return m;
}

constexpr Mat4x4f identity_Mat4x4f()
{
    return identity<Mat4x4f>();
}

constexpr Mat4x4d identity_Mat4x4d()
{
    return identity<Mat4x4d>();
}

constexpr Mat3x3f identity_Mat3x3f()
{
    return identity<Mat3x3f>();
}

constexpr Mat3x3d identity_Mat3x3d()
{
    return identity<Mat3x3d>();
}

/// Create a rotation matrix. Rotates about one of the major axes.
Mat4x4d rotation_Mat4x4d(Axis axis, float angle);
/// Create a translation matrix
Mat4x4d translation_Mat4x4d(const Vec3d&);
/// Create a scaling matrix.
Mat4x4d scaling_Mat4x4d(const Vec3d&);

/// Create a rotation matrix. Rotates about one of the major axes.
Mat4x4f rotation_Mat4x4f(Axis axis, float angle);
/// Create a translation matrix
Mat4x4f translation_Mat4x4f(const Vec3f& v);
/// Create a scaling matrix.
Mat4x4f scaling_Mat4x4f(const Vec3f& v);

/// Create a rotation matrix. Rotates about one of the major axes.
Mat3x3d rotation_Mat3x3d(Axis axis, double angle);
/// Create a scaling matrix.
Mat3x3d scaling_Mat3x3d(const Vec3d& v);

/// Create a rotation matrix. Rotates about one of the major axes.
Mat3x3f rotation_Mat3x3f(Axis axis, float angle);
/// Create a scaling matrix.
Mat3x3f scaling_Mat3x3f(const Vec3f& v);

/// Compute inverse assuming that the upper-left 3x3 sub-matrix is
/// orthonormal (which is the case if the transformation is only
/// a concatenation of rotations and translations).
constexpr Mat4x4d invert_ortho(const Mat4x4d& m)
{
    Vec3d rx(m[0][0], m[1][0], m[2][0]);
    Vec3d ry(m[0][1], m[1][1], m[2][1]);
    Vec3d rz(m[0][2], m[1][2], m[2][2]);
    Vec3d t(m[0][3], m[1][3], m[2][3]);

    return Mat4x4d(Vec4d(rx, -dot(t, rx)),
                   Vec4d(ry, -dot(t, ry)),
                   Vec4d(rz, -dot(t, rz)),
                   Vec4d(0.0, 0.0, 0.0, 1.0));
}

/// Creates a perspective projection similar to gluPerspective
/// Description from gluPerspective: perspective_Mat4x4f specifies a viewing frustum into
/// the world coordinate system. In general, the aspect ratio in perspective_Mat4x4f
/// should match the aspect ratio of the associated viewport. For example, aspect = 2.0
/// means the viewer's angle of view is twice as wide in x as it is in y. If the viewport
/// is twice as wide as it is tall, it displays the image without distortion.
Mat4x4f perspective_Mat4x4f(float fovy, float aspect, float zNear, float zFar);

/// Creates a perspective matrix similar to glFrustum
Mat4x4f frustum_Mat4x4f(float left,
                        float right,
                        float bottom,
                        float top,
                        float nearVal,
                        float farVal);

/// Creates an orthographic projection matrix (similar to glOrtho)
Mat4x4f ortho_Mat4x4f(float left,
                      float right,
                      float bottom,
                      float top,
                      float nearVal,
                      float farVal);

/// Creates a 2D orthographic projection matrix (similar to gluOrtho2D)
Mat4x4f ortho2D_Mat4x4f(float left, float right, float bottom, float top);

/// Creates a view matrix similar to gluLookAt
Mat4x4f lookAt_Mat4x4f(const Vec3f& eye, const Vec3f& at, const Vec3f& up);

/** Compute inverse assuming that the upper-left 3x3 sub-matrix is
    orthonormal (which is the case if the transformation is only
    a concatenation of rotations and translations).
*/
constexpr Mat4x4f invert_ortho(const Mat4x4f& m)
{
    Vec3f rx(m[0][0], m[1][0], m[2][0]);
    Vec3f ry(m[0][1], m[1][1], m[2][1]);
    Vec3f rz(m[0][2], m[1][2], m[2][2]);
    Vec3f t(m[0][3], m[1][3], m[2][3]);

    return Mat4x4f(Vec4f(rx, -dot(t, rx)),
                   Vec4f(ry, -dot(t, ry)),
                   Vec4f(rz, -dot(t, rz)),
                   Vec4f(0.0f, 0.0f, 0.0f, 1.0f));
}
}

#endif //GEL_MAT_H
