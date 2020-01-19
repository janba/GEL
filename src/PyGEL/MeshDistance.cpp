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

class MeshDistance {
    Geometry::AABBTree aabb_tree;
public:
    MeshDistance(HMesh::Manifold* m);
    
    float signed_distance(const CGLA::Vec3f& p, float upper);
    bool ray_inside_test(const CGLA::Vec3f& p, int no_rays);
};


MeshDistance::MeshDistance(Manifold* m) {
    build_AABBTree(*m, aabb_tree);
}

float MeshDistance::signed_distance(const CGLA::Vec3f& p, float upper){
    
    return aabb_tree.compute_signed_distance(p);
    
}

bool MeshDistance::ray_inside_test(const CGLA::Vec3f& p, int no_rays) {
    auto rand_vec = []() {return Vec3f(gel_rand()/double(GEL_RAND_MAX),
                                       gel_rand()/double(GEL_RAND_MAX),
                                       gel_rand()/double(GEL_RAND_MAX));
    };
    
    int even=0;
    int odd=0;
    for (int i=0;i<no_rays;++i) {
        int cnt = aabb_tree.intersect_cnt(p, rand_vec());
        if(cnt % 2 == 0)
            ++even;
        else
            ++odd;
    }
    return odd > even;
}


// -----------------------------------------
// C API
// -----------------------------------------

MeshDistance_ptr MeshDistance_new(Manifold_ptr _m) {
    Manifold* m = reinterpret_cast<Manifold*>(_m);
    return reinterpret_cast<MeshDistance_ptr>(new MeshDistance(m));
}

void MeshDistance_delete(MeshDistance_ptr _self) {
    MeshDistance* self = reinterpret_cast<MeshDistance*>(_self);
    delete self;
}

float MeshDistance_signed_distance(MeshDistance_ptr _self,
                      const float* _p,
                      float upper) {
    MeshDistance* self = reinterpret_cast<MeshDistance*>(_self);
    const Vec3f* p = reinterpret_cast<const Vec3f*>(_p);
    return self->signed_distance(*p, upper);
}

bool MeshDistance_ray_inside_test(MeshDistance_ptr _self,
                                  const float* _p,
                                  int no_rays) {
    MeshDistance* self = reinterpret_cast<MeshDistance*>(_self);
    const Vec3f* p = reinterpret_cast<const Vec3f*>(_p);
    return self->ray_inside_test(*p, no_rays);
}

