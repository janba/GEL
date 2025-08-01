/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithQuat.h
 * Abstract quaternion class
 */

#ifndef CGLA_ARITHQUAT_H
#define CGLA_ARITHQUAT_H

#include <GEL/CGLA/CGLA-util.h>
#include <GEL/CGLA/ArithVecFloat.h>
#include <GEL/CGLA/ArithSqMatFloat.h>
#include <cmath>

namespace CGLA
{
/** \brief A generic quaternion class.

  Quaternions are algebraic entities useful for rotation. */

template <class T, class V, class Q>
class ArithQuat {
public:
    /// Vector part of quaternion
    V qv;

    /// Scalar part of quaternion
    T qw;

    /// Default construct quaternion
    constexpr ArithQuat() = default;

    /// Construct quaternion from vector and scalar
    constexpr explicit ArithQuat(const V& imaginary, T real = 1.0f) : qv(imaginary), qw(real) {}

    /// Construct quaternion from four scalars
    constexpr ArithQuat(T x, T y, T z, T _qw) : qv(x, y, z), qw(_qw) {}

    /// @name Accessors
    /// @{

    /// Assign values to a quaternion
    constexpr void set(const V& imaginary, T real = 1.0f)
    {
        qv = imaginary;
        qw = real;
    }

    /// Assign values to a quaternion
    constexpr void set(T x, T y, T z, T _qw)
    {
        qv.set(x, y, z);
        qw = _qw;
    }

    /*constexpr void set(const Vec4f& v)
    {
      qv.set(v[0], v[1], v[2]);
      qw = v[3];
    }*/

    /// Get values from a quaternion
    constexpr void get(T& x, T& y, T& z, T& _qw) const
    {
        x = qv[0];
        y = qv[1];
        z = qv[2];
        _qw = qw;
    }

    /// Get imaginary part of a quaternion
    constexpr V get_imaginary_part() const { return qv; }

    /// Get real part of a quaternion
    constexpr T get_real_part() const { return qw; }

    /// @}
    /// @name Rotation operations
    /// @{

    /// Obtain angle of rotation and axis
    constexpr void get_rot(T& angle, V& v)
    {
        angle = 2 * std::acos(qw);

        if (angle < TINY)
            v = V(static_cast<T>(1.0), static_cast<T>(0.0), static_cast<T>(0.0));
        else
            v = qv * (static_cast<T>(1.0) / std::sin(angle));

        if (angle > M_PI)
            v = -v;

        v.normalize();
    }

    /// Construct a Quaternion from an angle and axis of rotation.
    constexpr void make_rot(T angle, const V& v)
    {
        angle /= static_cast<T>(2.0);
        qv = CGLA::normalize(v) * std::sin(angle);
        qw = std::cos(angle);
    }

    /** Construct a Quaternion rotating from the direction given
	by the first argument to the direction given by the second.*/
    constexpr void make_rot(const V& s, const V& t)
    {
        T tmp = std::sqrt(2 * (1 + dot(s, t)));
        qv = cross(s, t) * (1.0 / tmp);
        qw = tmp / 2.0;
    }

    /// Construct a Quaternion from a rotation matrix.
    template <class VT, class MT, unsigned int ROWS>
    constexpr void make_rot(ArithSqMatFloat<VT, MT, ROWS>& m)
    {
        static_assert(ROWS == 3 || ROWS == 4);
        T trace = m[0][0] + m[1][1] + m[2][2] + static_cast<T>(1.0);

        //If the trace of the matrix is greater than zero, then
        //perform an "instant" calculation.
        if (trace > TINY) {
            T S = sqrt(trace) * static_cast<T>(2.0);
            qv = V(m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1]);
            qv /= S;
            qw = static_cast<T>(0.25) * S;
        } else {
            //If the trace of the matrix is equal to zero (or negative...) then identify
            //which major diagonal element has the greatest value.
            //Depending on this, calculate the following:

            if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
                // Column 0:
                T S = sqrt(static_cast<T>(1.0) + m[0][0] - m[1][1] - m[2][2]) * static_cast<T>(2.0);
                qv[0] = 0.25f * S;
                qv[1] = (m[1][0] + m[0][1]) / S;
                qv[2] = (m[0][2] + m[2][0]) / S;
                qw = (m[2][1] - m[1][3]) / S;
            } else if (m[1][1] > m[2][2]) {
                // Column 1:
                T S = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
                qv[0] = (m[1][0] + m[0][1]) / S;
                qv[1] = 0.25 * S;
                qv[2] = (m[2][1] + m[1][2]) / S;
                qw = (m[0][2] - m[2][0]) / S;
            } else {
                // Column 2:
                T S = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
                qv[0] = (m[0][2] + m[2][0]) / S;
                qv[1] = (m[2][1] + m[1][2]) / S;
                qv[2] = 0.25 * S;
                qw = (m[1][0] - m[0][1]) / S;
            }
        }
        //The quaternion is then defined as:
        //  Q = | X Y Z W |
    }

    /// @}
    /// @name Binary operators
    /// @{

    /// Partial equality
    constexpr bool operator==(const ArithQuat<T, V, Q>& q) const
    {
        return qw == q.qw && qv == q.qv;
    }

    /// Partial inequality
    constexpr bool operator!=(const ArithQuat<T, V, Q>& q) const
    {
        return qw != q.qw || qv != q.qv;
    }

    /// Multiply two quaternions. (Combine their rotation)
    constexpr Q operator*(const ArithQuat<T, V, Q>& q) const
    {
        return Q(cross(qv, q.qv) + qv * q.qw + q.qv * qw,
                 qw * q.qw - dot(qv, q.qv));
    }

    /// Multiply scalar onto quaternion.
    constexpr Q operator*(T scalar) const
    {
        return Q(qv * scalar, qw * scalar);
    }

    /// Add two quaternions.
    constexpr Q operator+(const ArithQuat<T, V, Q>& q) const
    {
        return Q(qv + q.qv, qw + q.qw);
    }

    /// @}
    /// @name Unary operators
    /// @{

    /// Compute the additive inverse of the quaternion
    constexpr Q operator-() const { return Q(-qv, -qw); }

    /// Compute norm of quaternion
    constexpr T norm() const { return dot(qv, qv) + qw * qw; }

    /// Return conjugate quaternion
    constexpr Q conjugate() const { return Q(-qv, qw); }

    /// Compute the multiplicative inverse of the quaternion
    constexpr Q inverse() const { return Q(conjugate() * (1 / norm())); }

    /// Normalize quaternion.
    constexpr Q normalize() { return Q((*this) * (1 / norm())); }

    /// @}
    /// @name Application
    /// @{

    /// Rotate vector according to quaternion
    constexpr V apply(const V& vec) const
    {
        return ((*this) * Q(vec) * inverse()).qv;
    }

    /// Rotate vector according to unit quaternion
    constexpr V apply_unit(const V& vec) const
    {
        assert(fabs(norm() - 1.0) < SMALL);
        return ((*this) * Q(vec) * conjugate()).qv;
    }
};

/// @name Non-member ops
/// @related ArithQuat
/// @{

/// @brief Multiply with scalar
template <class T, class V, class Q>
constexpr Q operator*(T scalar, const ArithQuat<T, V, Q>& q)
{
    return q * scalar;
}

/** Perform linear interpolation of two quaternions.
    The last argument is the parameter used to interpolate
    between the two first. SLERP - invented by Shoemake -
    is a good way to interpolate because the interpolation
    is performed on the unit sphere.
*/
template <class T, class V, class Q>
constexpr Q slerp(const ArithQuat<T, V, Q>& q0, const ArithQuat<T, V, Q>& q1, T t)
{
    T angle = std::acos(q0.qv[0] * q1.qv[0] + q0.qv[1] * q1.qv[1]
        + q0.qv[2] * q1.qv[2] + q0.qw * q1.qw);
    return (q0 * std::sin((1 - t) * angle) + q1 * std::sin(t * angle)) * (1 / std::sin(angle));
}

/// @}
/// @name stdio operations
/// @related ArithQuat
/// @{

/// Print quaternion to stream.
template <class T, class V, class Q>
std::ostream& operator<<(std::ostream& os, const ArithQuat<T, V, Q>& v)
{
    os << "[ ";
    for (unsigned int i = 0; i < 3; i++) os << v.qv[i] << " ";
    os << "~ " << v.qw << " ";
    os << "]";

    return os;
}

///@}

}


#endif // CGLA_ARITHQUAT_H
