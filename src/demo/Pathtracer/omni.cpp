#include "omni.h"

using namespace CGLA;

#include "scene.h"

extern scene* current;

omni::omni(const CGLA::Vec3f& phi)
{
    set_power(phi);
}

void omni::set_power(const CGLA::Vec3f& phi)
{
    power_ = phi;
}

bool omni::sample(const ray& r,
                  const hit_info& hi, 
                  Vec3f& Li,
                  Vec3f& w) const
{
    Vec3f light_point(trs_[0][3], trs_[1][3], trs_[2][3]);

    w = normalize(light_point - hi.position);
    float r2 = sqr_length(light_point - hi.position);

    //send shadow feeler
    ray s;
    s.origin = hi.position + epsilon * w;
    s.direction = w;
    s.distance = std::sqrt(r2) - 2.f * epsilon;
    s.depth = r.depth + 1;

    if (current->intersect(s))
        return false;

    //compute radiance
    Li = power_ / (4.f * float(M_PI) * r2);

    return true;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
