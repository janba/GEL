#include "metal.h"

using namespace CGLA;

metal::metal(const CGLA::Vec3f& c, float s, float ior, float ext)
: material(ior, ext), color_(c), shininess_(s)
{
}

void metal::sample(const ray& r, hit_info& hi) const
{
    //compute fresnel
    float cosi = dot(-r.direction, hi.shading_normal);

    if (cosi <= 0.f)
        return;

    float f = fresnel_conductor(cosi, ior_, extinction_);

    if (f < 0.f || f > 1.f)
        return;

    if (shininess_ < 1000.f)
    {
        hi.glossy = f * color_;
        hi.shininess = shininess_;
    }
    else
    {
        hi.reflection = f * color_;
    }

    hi.ior = ior_;
    hi.extinction = extinction_;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
