#ifndef PATHTRACER_GLASS_H
#define PATHTRACER_GLASS_H

#include "material.h"

class glass : public material
{
public:
    glass(const CGLA::Vec3f& c, float ior);

    void sample(const ray&, hit_info&) const;

private:
    CGLA::Vec3f color_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
