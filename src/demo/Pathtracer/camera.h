#ifndef __PATHTRACER_CAMERA__H__
#define __PATHTRACER_CAMERA__H__

#pragma once

#include <GEL/CGLA/Vec2f.h>

#include "core.h"

class camera
{
public:
    camera(const CGLA::Vec3f& eye, const CGLA::Vec3f& center,
           const CGLA::Vec3f& up, float focal);

    ray generate(const CGLA::Vec2f&) const;

private:
    CGLA::Vec3f eye_;
    CGLA::Vec3f x_, y_, z_;
    float focal_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
