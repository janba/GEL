//
//  MeshDistance.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 13/12/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef MeshDistance_hpp
#define MeshDistance_hpp

#include <vector>
#include <utility>
#include <GEL/Geometry/build_bbtree.h>
#include "Manifold.h"

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

    MeshDistance* MeshDistance_new(Manifold* m);
    void MeshDistance_delete(MeshDistance* self);

    std::vector<float> MeshDistance_signed_distance(MeshDistance* self, const std::vector<float>& p, float upper);
    std::vector<int> MeshDistance_ray_inside_test(MeshDistance* self, const std::vector<float>& p, int no_rays);
    bool MeshDistance_ray_intersect(MeshDistance* self, CGLA::Vec3f& _p, CGLA::Vec3f& _d, float *t);
}

#endif /* MeshDistance_hpp */