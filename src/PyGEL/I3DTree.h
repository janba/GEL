//
//  I3DTree.h
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 02/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef I3DTree_h
#define I3DTree_h

#include <vector>
#include <utility>
#include "IntVector.h"
#include "Vec3dVector.h"
#include <GEL/CGLA/Vec3d.h>
#include <GEL/Geometry/KDTree.h>

// Define I3DTree as a specific KDTree type
using I3DTree = Geometry::KDTree<CGLA::Vec3d, size_t>;

namespace PyGEL {
    using I3DTreePtr = I3DTree*;
    using I3DTree_ptr = I3DTreePtr; // C-style alias
    
    I3DTree_ptr I3DTree_new();
    void I3DTree_delete(I3DTree_ptr self);
    void I3DTree_insert(I3DTree_ptr tree, double x, double y, double z, size_t v);
    void I3DTree_build(I3DTree_ptr tree);
    std::pair<std::vector<double>, size_t> I3DTree_closest_point(I3DTree_ptr tree, double x, double y, double z, double r);
    std::pair<Vec3dVector, IntVector> I3DTree_in_sphere(I3DTree_ptr tree, double x, double y, double z, double r);
    std::pair<Vec3dVector, IntVector> I3DTree_m_closest_points(I3DTree_ptr tree, double x, double y, double z, double r, int m);
}

#endif /* I3DTree_h */
