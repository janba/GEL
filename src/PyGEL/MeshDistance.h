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
#include "Manifold.h"

// Forward declaration for MeshDistance - this is a PyGEL-specific class
class MeshDistance;

namespace PyGEL {
    using MeshDistancePtr = MeshDistance*;
    using MeshDistance_ptr = MeshDistancePtr; // C-style alias
    
    MeshDistance_ptr MeshDistance_new(Manifold_ptr m);
    void MeshDistance_delete(MeshDistance_ptr self);
    
    std::vector<float> MeshDistance_signed_distance(MeshDistance_ptr self, const std::vector<float>& p, float upper);
    std::vector<int> MeshDistance_ray_inside_test(MeshDistance_ptr self, const std::vector<float>& p, int no_rays);
    std::pair<bool, float> MeshDistance_ray_intersect(MeshDistance_ptr self, const std::vector<float>& p, const std::vector<float>& d);
}

#endif /* MeshDistance_hpp */