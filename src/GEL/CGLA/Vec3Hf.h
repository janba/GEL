/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Vec3Hf.h
 * @brief Homogeneous 3D vector class - really a 4D vector.
 */

#ifndef __CGLA_VEC3HF_H__
#define __CGLA_VEC3HF_H__

#include "Vec4f.h"

namespace CGLA {

    /** A 3D homogeneous vector is simply a four D vector. */

    class Vec3Hf: public Vec4f {

    public:
        Vec3Hf();

        Vec3Hf(float _a);

        Vec3Hf(float _a, float _b, float _c, float _d);

        Vec3Hf(float _a, float _b, float _c);

        Vec3Hf(const Vec3f &v);

        Vec3Hf(const Vec3f &v, float _d);

        /// Divide all coordinates by the fourth coordinate
        inline void de_homogenize()
        {
            float k = 1.0f/(*this)[3];
            (*this)[0] = (*this)[0]*k;
            (*this)[1] = (*this)[1]*k;
            (*this)[2] = (*this)[2]*k;
            (*this)[3] = 1.0f;
        }
    };

    Vec3Hf::Vec3Hf() : Vec4f() {}

    Vec3Hf::Vec3Hf(float _a) : Vec4f(_a) {}

    Vec3Hf::Vec3Hf(float _a, float _b, float _c, float _d) : Vec4f(_a, _b, _c, _d) {}

    Vec3Hf::Vec3Hf(float _a, float _b, float _c) : Vec4f(_a, _b, _c) {}

    Vec3Hf::Vec3Hf(const Vec3f &v) : Vec4f(v) {}

    Vec3Hf::Vec3Hf(const Vec3f &v, float _d) : Vec4f(v, _d) {}

}
#endif

