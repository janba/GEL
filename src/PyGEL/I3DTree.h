//
//  PyGEL.h
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 02/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef PyGEL_h
#define PyGEL_h

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <GEL/CGLA/Vec3d.h>
#include <GEL/Geometry/KDTree.h>

#include "IntVector.h"
#include "Vec3dVector.h"

template class Geometry::KDTree<CGLA::Vec3d, size_t>;
using I3DTree = Geometry::KDTree<CGLA::Vec3d, size_t>;

extern "C" {
    DLLEXPORT I3DTree* I3DTree_new();
    DLLEXPORT void I3DTree_delete(I3DTree* self);
    DLLEXPORT void I3DTree_insert(I3DTree* tree, double x, double y, double z, size_t v);
    DLLEXPORT void I3DTree_build(I3DTree* tree);
    DLLEXPORT size_t I3DTree_closest_point(I3DTree* tree, double x, double y, double z, double r,
                                    CGLA::Vec3d* key, size_t* val);
    DLLEXPORT size_t I3DTree_in_sphere(I3DTree* tree, double x, double y, double z, double r,
                             Vec3dVector* keys, IntVector* vals);
}

#endif /* PyGEL_h */
