/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3usi.h
 * @brief 3D unsigned short integer vector class.
 */

#ifndef CGLA_VEC3USI_H
#define CGLA_VEC3USI_H

#include <GEL/CGLA/Vec.h>

namespace CGLA
{


// /** \brief Unsigned short int 3D vector class.
//
//         This class is mainly useful if we need a 3D int vector that takes up
//         less room than a Vec3i but holds larger numbers than a Vec3c. */
// class Vec3usi : public ArithVec3Int<int, Vec3usi> {
// public:
//     /// Construct 0 vector.
//     constexpr Vec3usi() = default;
//
//     /// Construct a Vec3usi
//     constexpr Vec3usi(USInt _a, USInt _b, USInt _c):
//         ArithVec3Int(_a, _b, _c) {}
//
//     /// Construct a Vec3usi from a Vec3i.
//     constexpr explicit Vec3usi(const Vec3i& v):
//         ArithVec3Int(v[0] & 0xffff, v[1] & 0xffff, v[2] & 0xffff) {}
// };
}
#endif
