#ifndef __PATHTRACER_LUMINAIRE__H__
#define __PATHTRACER_LUMINAIRE__H__

#include <GEL/CGLA/Mat4x4f.h>

#include "core.h"

class luminaire
{
public:
    luminaire(int samples = 1);
    virtual ~luminaire(void) = 0;

    //transformation matrix
    const CGLA::Mat4x4f& transform(void) const;
    void set_transform(const CGLA::Mat4x4f&);

    //sampling
    virtual bool sample(const ray&, const hit_info&, CGLA::Vec3f& L,
                        CGLA::Vec3f& w) const = 0;

    //sample count for direct lighting integration
    int samples(void) const;
    void set_samples(int);

protected:
    int samples_;
    CGLA::Mat4x4f trs_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
