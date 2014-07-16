#ifndef __PATHTRACER_CORE__H__
#define __PATHTRACER_CORE__H__

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3f.h>

#include <limits>

struct ray
{
    ray(void) : origin(0.f), direction(0.f), distance(0.f), depth(0) {}

    CGLA::Vec3f origin;
    CGLA::Vec3f direction;

    float distance;
    int depth;
};

struct hit_info
{
    hit_info(void) :
        distance(0.f),
        position(0.f), geometric_normal(0.f), shading_normal(0.f),
        inside(false),
        emitted(0.f),
        diffuse(0.f), glossy(0.f), shininess(0.f),
        reflection(0.f), refraction(0.f),
        ior(0.f), extinction(0.f) {}

    float distance;
    CGLA::Vec3f position;
    CGLA::Vec3f geometric_normal;
    CGLA::Vec3f shading_normal;

    //The tangent, bitangent, and shading_normal form an orthonormal basis
    CGLA::Vec3f tangent;
    CGLA::Vec3f bitangent;

    //texture coordinates
    CGLA::Vec2f texcoords;

    //flag to indicate if the inside of the object (wrt the normal) was hit
    bool inside;

    //emitted radiance
    CGLA::Vec3f emitted;

    //diffuse reflectivity (rho)
    CGLA::Vec3f diffuse;

    //parameters for the modified Phong reflectance model
    CGLA::Vec3f glossy;
    float shininess;

    //parameters for perfect specular reflection/refraction
    CGLA::Vec3f reflection;
    CGLA::Vec3f refraction;

    //index of refraction + complex part
    float ior;
    float extinction;
};

//RGB to intensity
inline float intensity(const CGLA::Vec3f& v)
{
    return v[0] * 0.27f + v[1] * 0.67f + v[2] * 0.06f;
}

//standard fresnel formula (see eg. Jensen p. 23)
inline float fresnel_dielectric(float cosi, float cost, float eta)
{
    assert(cosi>0.f && cost>0.f);

    float etai = 1.f;
    float etat = eta;

    float Fs = (etai*cosi - etat*cost) / (etai*cosi + etat*cost);
    float Fp = (etat*cosi - etai*cost) / (etat*cosi + etai*cost);

    float F = 0.5f * (Fs*Fs + Fp*Fp);
    assert(F>=0.f && F<=1.f);
    return F;
}

//standard fresnel formula (see eg. Pharr et all p. 422)
inline float fresnel_conductor(float cosi, float eta, float extinction)
{
    assert(cosi >= 0.f);

    float z = eta*eta + extinction*extinction;

    float Fp = (z * cosi*cosi - 2.f * eta * cosi + 1.f) /
        (z + 2.f * eta * cosi + 1.f);

    float Fs = (z - 2.f * eta * cosi + cosi * cosi) /
        (z + 2.f * eta * cosi + cosi * cosi);

    float F = (Fp + Fs) * 0.5f;
    assert(F>=0.f && F<=1.f);
    return F;
}

//reflect vector around normal. n should be normalized.
inline CGLA::Vec3f reflect(const CGLA::Vec3f& n, const CGLA::Vec3f& r)
{
    return 2.f * n * dot(n, r) - r;
}

//eta = eta_from / eta_to
inline bool refract(const CGLA::Vec3f & n, const CGLA::Vec3f & i, 
                    float eta, CGLA::Vec3f & t)
{
    float c1 = dot(i,n);
    float c2 = 1.f - eta*eta * (1.f - c1*c1);

    if (c2 < 0.f)
        return false;

    float c3 = std::sqrt(c2);
    t = -eta*i + (eta*c1 - c3) * n;

    return true;
}

//helper function, clamps v to range [i; a]
template <class T>
T clamp(const T& v, const T& i, const T& a)
{
    return std::min(std::max(v, i), a);
}

//mersenne twister
double genrand_real2(void);
inline float mt_random(void)
{
    return float(genrand_real2());
}

//evaluation of the diffuse brdf
inline CGLA::Vec3f lambertian_brdf(const hit_info& hi,
                                   const CGLA::Vec3f& wi,
                                   const CGLA::Vec3f& wo)
{
    return hi.diffuse / float(M_PI);
}

//sampling of diffuse brdf, wi = sampled direction, should return probability
inline float sample_lambertian(const hit_info& hi,
                            const CGLA::Vec3f& wo, CGLA::Vec3f& wi)
{
    CGLA::Vec3f x, y, z = hi.shading_normal;
    orthogonal(z, x, y);
    float e0 = mt_random();
    float e1 = mt_random();
    float cost = std::sqrt(e0);
    float sint = std::sqrt(1.f - e0);
    float cosp = std::cos(e1 * 2.f * float(M_PI));
    float sinp = std::sin(e1 * 2.f * float(M_PI));
    CGLA::Vec3f dir(sint*cosp, sint*sinp, cost);

    wi = dir(0) * x + dir(1) * y + dir(2) * z;
    return 1.f / float(M_PI) * cost;
}

//evaluation of the phong brdf
inline CGLA::Vec3f phong_brdf(const hit_info& hi,
                               const CGLA::Vec3f& wi,const CGLA::Vec3f& wo)
{
    float dwi = dot(wi, hi.geometric_normal);
    float dwo = dot(wo, hi.geometric_normal);
    float same_hemisphere = dwi * dwo > 0.f;

    if (same_hemisphere)
    {
        //glossy contribution
        if (intensity(hi.glossy) > 0.f)
        {
            float cost = std::max(dot(reflect(hi.shading_normal, wi), wo), 0.f);
            float cosn = std::pow(cost, hi.shininess);
            return hi.glossy * (hi.shininess+2.f)/(2.f * float(M_PI))*cosn;
        }
    }
    return CGLA::Vec3f(0.f);
}

//sampling of phong brdf, wi = sampled direction, should return probability
inline float sample_phong(const hit_info& hi,
                          const CGLA::Vec3f& wo, CGLA::Vec3f& wi)
{
    float e0 = mt_random();
    float e1 = mt_random();

    float cost = std::pow(e0, 1.f/(hi.shininess + 1.f));
    float sint = std::sqrt(1.f - std::pow(e0, 2.f/(hi.shininess + 1.f)));
    float cosp = std::cos(e1 * 2.f * float(M_PI));
    float sinp = std::sin(e1 * 2.f * float(M_PI));

    CGLA::Vec3f dir(sint*cosp, sint*sinp, cost);

    CGLA::Vec3f x, y, z = reflect(hi.shading_normal, wo);
    orthogonal(z, x, y);

    float cosn = std::pow(dir(2), hi.shininess);
    wi = dir(0) * x + dir(1) * y + dir(2) * z;

    return (hi.shininess+1.f) / (2.f * float(M_PI)) * cosn;
}

//convenience function to evaluate the non-specular parts of the bsdf
inline CGLA::Vec3f bsdf_evaluate(const hit_info& hi,
                                 const CGLA::Vec3f& wi, const CGLA::Vec3f& wo)
{
    CGLA::Vec3f fs(0.f);

    fs += lambertian_brdf(hi, wi, wo);
    fs += phong_brdf(hi, wi, wo);

    return fs;
}

//epsilon for shadow testing
static const float epsilon = 1e-5f;

static const float step = 1e-1f;

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
