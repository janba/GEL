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

#include <stdbool.h>

#include "IntVector.h"
#include "Vec3dVector.h"

typedef char* I3DTree_ptr;

#ifdef __cplusplus
extern "C" {
#endif
    DLLEXPORT I3DTree_ptr I3DTree_new();
    DLLEXPORT void I3DTree_delete(I3DTree_ptr self);
    DLLEXPORT void I3DTree_insert(I3DTree_ptr tree, double x, double y, double z, size_t v);
    DLLEXPORT void I3DTree_build(I3DTree_ptr tree);
    DLLEXPORT size_t I3DTree_closest_point(I3DTree_ptr tree, double x, double y, double z, double r,
                                           double* key, size_t* val);
    DLLEXPORT size_t I3DTree_in_sphere(I3DTree_ptr tree, double x, double y, double z, double r,
                                       Vec3dVector_ptr keys, IntVector_ptr vals);
    DLLEXPORT size_t I3DTree_m_closest_points(I3DTree_ptr tree, double x, double y, double z, double r, int m,
                                              Vec3dVector_ptr keys, IntVector_ptr vals);
    
#ifdef __cplusplus
}
#endif

#endif /* PyGEL_h */
