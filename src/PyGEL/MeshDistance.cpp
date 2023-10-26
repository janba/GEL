//
//  MeshDistance.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 13/12/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "MeshDistance.h"

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/build_bbtree.h>

using namespace HMesh;
using namespace CGLA;
using namespace Geometry;

class MeshDistance
{
    Geometry::AABBTree aabb_tree;

public:
    MeshDistance(HMesh::Manifold *m);

    float signed_distance(const CGLA::Vec3f &p, float upper);
    bool ray_inside_test(const CGLA::Vec3f &p, int no_rays);

    bool ray_intersect(CGLA::Vec3f &p, CGLA::Vec3f &d, float &t);
};

MeshDistance::MeshDistance(Manifold *m)
{
    build_AABBTree(*m, aabb_tree);
}

float MeshDistance::signed_distance(const CGLA::Vec3f &p, float upper)
{

    return aabb_tree.compute_signed_distance(p);
}

bool MeshDistance::ray_intersect(CGLA::Vec3f &p, CGLA::Vec3f &d, float &t)
{
    Ray ray(p, d);
    aabb_tree.intersect(ray);
    if (ray.has_hit)
    {
        t = ray.dist;
        p = ray.hit_pos;
        d = ray.hit_normal;
        return true;
    }
    return false;
}

bool MeshDistance::ray_inside_test(const CGLA::Vec3f &p, int no_rays)
{
    auto rand_vec = []()
    {
        return Vec3f(gel_rand() / double(GEL_RAND_MAX),
                     gel_rand() / double(GEL_RAND_MAX),
                     gel_rand() / double(GEL_RAND_MAX));
    };

    int even = 0;
    int odd = 0;
    for (int i = 0; i < no_rays; ++i)
    {
        int cnt = aabb_tree.intersect_cnt(p, rand_vec());
        if (cnt % 2 == 0)
            ++even;
        else
            ++odd;
    }
    return odd > even;
}

// -----------------------------------------
// C API
// -----------------------------------------

MeshDistance_ptr MeshDistance_new(Manifold_ptr _m)
{
    Manifold *m = reinterpret_cast<Manifold *>(_m);
    return reinterpret_cast<MeshDistance_ptr>(new MeshDistance(m));
}

void MeshDistance_delete(MeshDistance_ptr _self)
{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    delete self;
}

void MeshDistance_signed_distance(MeshDistance_ptr _self,
                                  int no_query_points,
                                  const float *p,
                                  float *d,
                                  float upper)
{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    for (int i = 0; i < no_query_points; ++i)
        d[i] = self->signed_distance(Vec3f(p[3 * i], p[3 * i + 1], p[3 * i + 2]), upper);
}

void MeshDistance_ray_inside_test(MeshDistance_ptr _self,
                                  int no_query_points,
                                  const float *p,
                                  int *s,
                                  int no_rays)
{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    for (int i = 0; i < no_query_points; ++i)
        s[i] = self->ray_inside_test(Vec3f(p[3 * i], p[3 * i + 1], p[3 * i + 2]), no_rays);
}

bool MeshDistance_ray_intersect(MeshDistance_ptr _self, float *_p, float *_d, float *t)
{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    Vec3f p(_p[0], _p[1], _p[2]);
    Vec3f d(_d[0], _d[1], _d[2]);
    return self->ray_intersect(p, d, *t);
    _p[0] = p[0];
    _p[1] = p[1];
    _p[2] = p[2];
    _d[0] = d[0];
    _d[1] = d[1];
    _d[2] = d[2];
}
