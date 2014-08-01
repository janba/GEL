#include "Vec3Hf.h"

namespace CGLA {
	Vec3Hf::Vec3Hf() : Vec4f() {}

	Vec3Hf::Vec3Hf(float _a) : Vec4f(_a) {}

	Vec3Hf::Vec3Hf(float _a, float _b, float _c, float _d) : Vec4f(_a, _b, _c, _d) {}

	Vec3Hf::Vec3Hf(float _a, float _b, float _c) : Vec4f(_a, _b, _c, 1.0f) {}

	Vec3Hf::Vec3Hf(const Vec3f &v) : Vec4f(v[0], v[1], v[2], 1.0f) {}

	Vec3Hf::Vec3Hf(const Vec3f &v, float _d) : Vec4f(v, _d) {}

	Vec3Hf::Vec3Hf(const Vec4f &v) : Vec4f(v) {}
}