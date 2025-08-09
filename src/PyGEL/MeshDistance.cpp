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

namespace PyGEL {

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

std::vector<float> MeshDistance_signed_distance(MeshDistance_ptr _self, const std::vector<float>& p, float upper)

{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    std::vector<float> d(p.size() / 3);
    for (size_t i = 0; i < p.size() / 3; ++i)
        d[i] = self->signed_distance(Vec3f(p[3 * i], p[3 * i + 1], p[3 * i + 2]), upper);
    return d;
}

std::vector<int> MeshDistance_ray_inside_test(MeshDistance_ptr self, const std::vector<float>& p, int no_rays)
{
    MeshDistance *dist = reinterpret_cast<MeshDistance *>(self);
    std::vector<int> results;
    for (size_t i = 0; i < p.size() / 3; ++i)
    {
        CGLA::Vec3f point(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
        results.push_back(dist->ray_inside_test(point, no_rays));
    }
    return results;
}


bool MeshDistance_ray_intersect(MeshDistance_ptr _self, CGLA::Vec3f& p, CGLA::Vec3f& d, float *t)
{
    MeshDistance *self = reinterpret_cast<MeshDistance *>(_self);
    return self->ray_intersect(p, d, *t);
}

} // namespace PyGEL
