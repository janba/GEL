#ifndef __PATHTRACER_MATTE__H__
#define __PATHTRACER_MATTE__H__

#include "material.h"

class matte : public material
{
public:
    //diffuse is rho (diffuse reflectance)
    matte(const CGLA::Vec3f& diffuse);

    void sample(const ray&, hit_info&) const;

private:
    CGLA::Vec3f diffuse_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
