#ifndef __PATHTRACER_OMNI__H__
#define __PATHTRACER_OMNI__H__

#include "luminaire.h"

class omni : public luminaire
{
public:
    omni(const CGLA::Vec3f& phi);

    void set_power(const CGLA::Vec3f& phi);

    bool sample(const ray&, const hit_info&, CGLA::Vec3f& L,
                CGLA::Vec3f& w) const;

private:

    //power in W
    CGLA::Vec3f power_;

};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
