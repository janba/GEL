#include "matte.h"

using namespace CGLA;

matte::matte(const CGLA::Vec3f& r) : material(1.f, 0.f), diffuse_(r)
{
    assert(diffuse_[0]>=0.f && diffuse_[0]<1.f);
    assert(diffuse_[1]>=0.f && diffuse_[1]<1.f);
    assert(diffuse_[2]>=0.f && diffuse_[2]<1.f);
}

void matte::sample(const ray& r, hit_info& hi) const
{
    hi.diffuse = diffuse_;
    hi.ior = hi.inside ? 1.f / ior_ : ior_;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
