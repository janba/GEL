//
// Created by Cem Akarsubasi on 8/4/25.
//

/// @file GEL/CGLA/Definitions.h
/// This is a definition file that replaces the 20+ separate classes in different headers that CGLA previously
/// included. Do not include this directly, instead include @ref GEL/CGLA/Vec.h and @ref GEL/CGLA/Mat.h

#ifndef CGLA_DEFINITIONS_H
#define CGLA_DEFINITIONS_H

#define GEL_STRINGIFY0(x) #x
#define GEL_STRINGIFY1(x) GEL_STRINGIFY0(x)
#ifdef __GNUC__
    #define GEL_COMPILER_WARNING_COMPOSE(x) GCC warning x
#else // Mainly MSVC
    #define GEL_COMPILER_WARNING_COMPOSE(x) message(x)
#endif
#define GEL_WARNING(x) _Pragma(GEL_STRINGIFY1(GEL_COMPILER_WARNING_COMPOSE(x)))

#define GEL_CGLA_VEC_DEPRECATED GEL_WARNING("This file is deprecated, include GEL/CGLA/Vec.h instead")

#include <concepts>
#include <cstdint>

namespace CGLA
{
template <std::floating_point Float>
class Vec2;
template <std::floating_point Float>
class Vec3;
template <std::floating_point Float>
class Vec4;

template <std::integral Integral>
class Vec2Integral;
template <std::integral Integral>
class Vec3Integral;
template <std::integral Integral>
class Vec4Integral;

template <std::floating_point Float>
class Mat2x2;
template <std::floating_point Float>
class Mat3x3;
template <std::floating_point Float>
class Mat4x4;
template <std::floating_point Float>
class Mat2x3;
template <std::floating_point Float>
class Mat3x2;

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;
using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;
using Vec4f = Vec4<float>;
using Vec4d = Vec4<double>;

using Vec2c = Vec2Integral<std::int8_t>;
using Vec2si = Vec2Integral<std::int16_t>;
using Vec2i = Vec2Integral<std::int32_t>;
using Vec2l = Vec2Integral<std::int64_t>;

using Vec2uc = Vec2Integral<std::uint8_t>;
using Vec2usi = Vec2Integral<std::uint16_t>;
using Vec2ui = Vec2Integral<std::uint32_t>;
using Vec2ul = Vec2Integral<std::uint64_t>;


using Vec3c = Vec3Integral<std::int8_t>;
/// @brief Unsigned short int 3D vector class.
///
///  This class is mainly useful if we need a 3D int vector that takes up
///  less room than a Vec3i but holds larger numbers than a Vec3c.
using Vec3si = Vec3Integral<std::int16_t>;
using Vec3i = Vec3Integral<std::int32_t>;
using Vec3l = Vec3Integral<std::int64_t>;

using Vec3uc = Vec3Integral<std::uint8_t>;
using Vec3usi = Vec3Integral<std::uint16_t>;
using Vec3ui = Vec3Integral<std::uint32_t>;
using Vec3ul = Vec3Integral<std::uint64_t>;

using Vec4c  = Vec4Integral<std::int8_t>;
using Vec4si = Vec4Integral<std::int16_t>;
using Vec4i  = Vec4Integral<std::int32_t>;
using Vec4l  = Vec4Integral<std::int64_t>;

using Vec4uc  = Vec4Integral<std::uint8_t>;
using Vec4usi = Vec4Integral<std::uint16_t>;
using Vec4ui  = Vec4Integral<std::uint32_t>;
using Vec4ul  = Vec4Integral<std::uint64_t>;

using Mat2x2f = Mat2x2<float>;
using Mat2x2d = Mat2x2<double>;

using Mat3x3f = Mat3x3<float>;
using Mat3x3d = Mat3x3<double>;

using Mat4x4f = Mat4x4<float>;
using Mat4x4d = Mat4x4<double>;

using Mat2x3f = Mat2x3<float>;
using Mat2x3d = Mat2x3<double>;

using Mat3x2f = Mat3x2<float>;
using Mat3x2d = Mat3x2<double>;

using USInt = std::uint16_t;
using UChar = std::uint8_t;
}

#endif // CGLA_DEFINITIONS_H
