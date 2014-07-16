#include "glass.h"

using namespace CGLA;

glass::glass(const CGLA::Vec3f& c, float ior)
: material(ior, 0.f), color_(c)
{
}

void glass::sample(const ray& r, hit_info& hi) const
{
    hi.ior = hi.inside ? 1.f / ior_ : ior_;

    //compute fresnel
    float f;
    float cosi = dot(-r.direction, hi.shading_normal);
    cosi = clamp(cosi, -1.f, 1.f);
    float sint = 1.f / hi.ior * std::sqrt(std::max(0.f, 1.f - cosi*cosi));

    if (sint > 1.f)
    {
        //total internal reflection
        f = 1.f;
    }
    else if (cosi < 0.)
    {
        f = 0.f;
    }
    else
    {
        float cost = std::sqrt(std::max(0.f, 1.f - sint * sint));
        cost = clamp(cost, -1.f, 1.f);
        f = fresnel_dielectric(cosi, cost, hi.ior);

        if (f < 0.f || f > 1.f)
            assert(false);
    }

    hi.reflection = f * color_;
    hi.refraction = (1.f - f) * color_;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
