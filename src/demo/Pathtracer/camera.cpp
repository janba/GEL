#include "camera.h"

using namespace CGLA;

camera::camera(const CGLA::Vec3f& eye,
               const CGLA::Vec3f& cnt,
               const CGLA::Vec3f& up,
               float focal)
{
    eye_ = eye;

    z_ = normalize(eye - cnt);

    x_ = normalize(cross(up, z_));
    y_ = cross(z_, x_);

    x_ *= 0.012f;
    y_ *= 0.012f;

    focal_ = focal;
}

ray camera::generate(const Vec2f& uv) const
{
    Vec2f xy = 2.f * uv - Vec2f(1.f);

    ray r;
    r.origin = eye_;
    r.direction = normalize(-focal_*z_ + xy[0]*x_ + xy[1]*y_);
    r.depth = 0;
    r.distance = std::numeric_limits<float>::infinity();
    return r;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
