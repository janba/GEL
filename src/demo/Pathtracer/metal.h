#ifndef __PATHTRACER_METAL__H__
#define __PATHTRACER_METAL__H__

#include "material.h"

class metal : public material
{
public:
    metal(const CGLA::Vec3f& color, float shininess, float ior, float ext);

    void sample(const ray&, hit_info&) const;

private:
    CGLA::Vec3f color_;
    float shininess_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
