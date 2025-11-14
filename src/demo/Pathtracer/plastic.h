#ifndef PATHTRACER_PLASTIC_H
#define PATHTRACER_PLASTIC_H

#include "material.h"

class plastic : public material
{
public:
    plastic(const CGLA::Vec3f& diffuse, const CGLA::Vec3f& glossy,
        float shininess);

    void sample(const ray&, hit_info&) const;

private:
    CGLA::Vec3f diffuse_;
    CGLA::Vec3f glossy_;
    float shininess_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
