#ifndef __PATHTRACER_X__H__
#define __PATHTRACER_X__H__

#include "scene.h"

class pathtracer
{
public:
    pathtracer(int width, int height, bool explicit_direct, int subsamples);

    CGLA::Vec3f compute_pixel(int w, int h);

    void set_scene(scene*);

private:

    CGLA::Vec3f trace(const ray& r, bool include_emitted);

    //resolution
    int width_;
    int height_;

    //the scene
    scene* scene_;

    //use explicit integration of direct illumination
    bool explicit_direct_;

    //misc
    int subsamples_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
