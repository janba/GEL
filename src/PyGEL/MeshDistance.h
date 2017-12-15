//
//  MeshDistance.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 13/12/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef MeshDistance_hpp
#define MeshDistance_hpp

#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/build_bbtree.h>

class MeshDistance {
    Geometry::AABBTree aabb_tree;
public:
    MeshDistance(HMesh::Manifold* m);
    
    float signed_distance(const CGLA::Vec3f& p, float upper);
};

extern "C" {
    MeshDistance* MeshDistance_new(HMesh::Manifold* m);
    void MeshDistance_delete(MeshDistance*);
    
    float MeshDistance_signed_distance(MeshDistance* self,
                          const CGLA::Vec3f* p,
                          float upper);
}

#endif /* MeshDistance_hpp */
