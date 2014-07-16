#include "plastic.h"

using namespace CGLA;

plastic::plastic(const CGLA::Vec3f& d, const CGLA::Vec3f& g, float s)
: material(1.6f, 0.f), diffuse_(d), glossy_(g), shininess_(s)
{
}

void plastic::sample(const ray& r, hit_info& hi) const
{
    float cosi = dot(-r.direction, hi.shading_normal);
    cosi = clamp(cosi, -1.f, 1.f);
    if (cosi <= 0.f)
        return;
    float sint = 1.f / ior_ * std::sqrt(std::max(0.f, 1.f - cosi*cosi));
    float cost = std::sqrt(std::max(0.f, 1.f - sint * sint));
    cost = clamp(cost, -1.f, 1.f);
    float f = fresnel_dielectric(cosi, cost, ior_);

    hi.diffuse = diffuse_;
    hi.glossy = f*glossy_;
    hi.shininess = shininess_;
    hi.ior = ior_;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
