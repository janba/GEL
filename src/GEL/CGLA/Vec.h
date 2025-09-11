//
// Created by Cem Akarsubasi on 8/4/25.
//

#ifndef CGLA_VEC_H
#define CGLA_VEC_H

#include <GEL/CGLA/Definitions.h>
#include <GEL/CGLA/ArithVec3Int.h>
#include <GEL/CGLA/ArithVec4Int.h>
#include <GEL/CGLA/ArithVec2Float.h>
#include <GEL/CGLA/ArithVec3Float.h>
#include <GEL/CGLA/ArithVec4Float.h>

namespace CGLA
{
template <std::floating_point Float>
class Vec2 final : public ArithVec2Float<Float, Vec2<Float>> {
public:
    constexpr Vec2() = default;

    constexpr Vec2(Float _a, Float _b):
        ArithVec2Float<Float, Vec2>(_a, _b) {}

    template <std::floating_point OtherFloat>
    constexpr explicit Vec2(const Vec2<OtherFloat>& v):
        ArithVec2Float<Float, Vec2>(v[0], v[1]) {}

    constexpr explicit Vec2(Float a):
        ArithVec2Float<Float, Vec2>(a, a) {}

    template <class T, class V, unsigned int N>
    constexpr explicit Vec2(const ArithVec<T, V, N>& v):
        ArithVec2Float<Float, Vec2>(static_cast<Float>(v[0]),
                                    static_cast<Float>(v[1])) {}
};

template <std::integral Integral>
class Vec2Integral final : public ArithVec<Integral, Vec2Integral<Integral>, 2> {
public:
    /// Construct 0 vector
    constexpr Vec2Integral() = default;

    /// Construct 2D int vector
    constexpr Vec2Integral(Integral _a, Integral _b): ArithVec<Integral, Vec2Integral, 2>(_a, _b) {}

    /// Construct a 2D integer vector with two identical coordinates.
    constexpr explicit Vec2Integral(Integral a): ArithVec<Integral, Vec2Integral, 2>(a, a) {}

    /// Convert from 2D float vector
    template <std::floating_point Float>
    constexpr explicit Vec2Integral(const Vec2<Float>& v) :
        ArithVec<Integral, Vec2Integral, 2>(static_cast<Integral>(floor(v[0])), static_cast<Integral>(floor(v[1]))) {}
};


/** \brief A 3 dimensional vector. */
template <std::floating_point Float>
class Vec3 final : public ArithVec3Float<Float, Vec3<Float>> {
public:
    /// Construct 0 vector
    constexpr Vec3() = default;

    /// Construct vector
    constexpr Vec3(Float a, Float b, Float c): ArithVec3Float<Float, Vec3>(a, b, c) {}

    /// Construct vector where all coords = a
    constexpr explicit Vec3(Float a):
        ArithVec3Float<Float, Vec3>(a, a, a) {}

    /// Convert from int vector
    template <std::integral Integral>
    constexpr explicit Vec3(const Vec3Integral<Integral>& v):
        ArithVec3Float<Float, Vec3>(v[0], v[1], v[2]) {}

    /// Convert from float vector
    template <std::floating_point OtherFloat>
    constexpr explicit Vec3(const Vec3<OtherFloat>& v):
        ArithVec3Float<Float, Vec3>(v[0], v[1], v[2]) {}

    /// Construct from a 4D float vector (skipping the last value)
    template <std::floating_point OtherFloat>
    constexpr explicit Vec3(const Vec4<OtherFloat>& v): ArithVec3Float<Float, Vec3>(v[0], v[1], v[2]) {}
};

/** \brief 3D integer vector.

    This class does not really extend the template
        and hence provides only the basic facilities of an ArithVec.
        The class is typically used for indices to 3D voxel grids. */
template <std::integral Integral>
class Vec3Integral final : public ArithVec3Int<Integral, Vec3Integral<Integral>> {
public:
    /// Construct 0 vector.
    constexpr Vec3Integral() = default;

    /// Construct a 3D integer vector.
    constexpr Vec3Integral(Integral _a, Integral _b, Integral _c): ArithVec3Int<Integral, Vec3Integral>(_a, _b, _c) {}

    /// Construct a 3D integer vector with 3 identical coordinates.
    constexpr explicit Vec3Integral(Integral a): ArithVec3Int<Integral, Vec3Integral>(a, a, a) {}

    /// Construct from a Vec3.
    template <std::floating_point Float>
    explicit Vec3Integral(const Vec3<Float>& v): ArithVec3Int<Integral, Vec3Integral>(Integral(std::floor(v[0])),
        Integral(std::floor(v[1])),
        Integral(std::floor(v[2]))) {}

    /// Construct from Vec3Integral with a smaller size.
    template <std::integral OtherIntegral> requires (sizeof(Integral) >= sizeof(OtherIntegral))
    explicit Vec3Integral(const Vec3Integral<OtherIntegral>& v) : ArithVec3Int<Integral, Vec3Integral>(v[0], v[1], v[2])
    {}
};

/** \brief A four dimensional floating point vector.

        This class is also used (via typedef) for
        homogeneous vectors.
*/
template <std::floating_point Float>
class Vec4 final : public ArithVec4Float<Float, Vec4<Float>> {
public:
    /// Construct a (0,0,0,0) homogenous Vector
    constexpr Vec4(): ArithVec4Float<Float, Vec4>(0.0, 0.0, 0.0, 0.0) {}

    /// Construct a (0,0,0,0) homogenous Vector
    constexpr explicit Vec4(Float _a): ArithVec4Float<Float, Vec4>(_a, _a, _a, _a) {}

    /// Construct a 4D vector
    constexpr Vec4(Float _a, Float _b, Float _c, Float _d):
        ArithVec4Float<Float, Vec4>(_a, _b, _c, _d) {}

    /// Construct a homogenous vector from a non-homogenous.
    constexpr explicit Vec4(const Vec3<Float>& v, Float _d):
        ArithVec4Float<Float, Vec4>(v[0], v[1], v[2], _d) {}
};

// 4D Integral vector.
//
//     This class does not really extend the template
//         and hence provides only the basic facilities of an ArithVec.
//         The class is typically used for indices to 4D voxel grids.
template <std::integral Integral>
class Vec4Integral final: public ArithVec4Int<Integral, Vec4Integral<Integral>> {
public:
    /// Construct 0 vector.
    constexpr Vec4Integral() = default;

    /// Construct a 4D integer vector.
    constexpr Vec4Integral(int _a, int _b, int _c, int _d): ArithVec4Int<int, Vec4Integral>(_a, _b, _c, _d) {}

    /// Construct a 4D integer vector with 4 identical coordinates.
    constexpr explicit Vec4Integral(int a): ArithVec4Int<int, Vec4Integral>(a, a, a, a) {}

    /// Construct from a floating point vector.
    template <std::floating_point Float>
    explicit Vec4Integral(const Vec4<Float>& v) : ArithVec4Int<Integral, Vec4Integral>(Integral(std::floor(v[0])),
        Integral(std::floor(v[1])),
        Integral(std::floor(v[2])), Integral(std::floor(v[3]))) {}

    /// Construct from an integral vector.
    template <std::integral OtherIntegral>
    explicit Vec4Integral(const Vec4Integral<OtherIntegral>& v) : ArithVec4Int<Integral, Vec4Integral>(
        v[0], v[1], v[2], v[3]) {}
};

template <typename VecOrMat, typename Predicate>
auto any(const VecOrMat& v, Predicate&& p) -> bool
{
    for (const auto& elem: v) {
        if (p(elem)) {
            return true;
        }
    }
    return false;
}

template <typename VecOrMat, typename Predicate>
auto all(const VecOrMat& v, Predicate&& p) -> bool
{
    for (const auto& elem: v) {
        if (!p(v)) {
            return false;
        }
    }
    return true;
}

template <>
struct VecT_to_MatT<Vec2d> {
    using MatT = Mat2x2d;
};

template <>
struct VecT_to_MatT<Vec4d> {
    using MatT = Mat4x4d;
};

template <>
struct VecT_to_MatT<Vec4f> {
    using MatT = Mat4x4f;
};

template <>
struct VecT_to_MatT<Vec3d> {
    using MatT = Mat3x3d;
};

template <>
struct VecT_to_MatT<Vec3f> {
    using MatT = Mat3x3f;
};
}


#endif // CGLA_VEC_H
