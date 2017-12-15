//
//  MeshDistance.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 13/12/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "MeshDistance.h"
using namespace HMesh;
using namespace Geometry;

MeshDistance::MeshDistance(Manifold* m) {
    build_AABBTree(*m, aabb_tree);
}

float MeshDistance::signed_distance(const CGLA::Vec3f& p, float upper){
    
    return aabb_tree.compute_signed_distance(p);
    
}

MeshDistance* MeshDistance_new(HMesh::Manifold* m) {
    return new MeshDistance(m);
}

void MeshDistance_delete(MeshDistance* self) {
    delete self;
}

float MeshDistance_signed_distance(MeshDistance* self,
                      const CGLA::Vec3f* p,
                      float upper) {
    return self->signed_distance(*p, upper);
}
