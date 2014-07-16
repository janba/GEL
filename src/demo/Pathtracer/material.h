#ifndef __PATHTRACER_MATERIAL__H__
#define __PATHTRACER_MATERIAL__H__

#include "core.h"

class material
{
public:
    material(float ior, float ext);
    virtual ~material(void) = 0;

    virtual void sample(const ray&, hit_info&) const = 0;

protected:

    //index of refraction (real part)
    float ior_;

    //index of refraction (complex part)
    float extinction_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
