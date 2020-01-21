//
//  Vec3dVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Vec3dVector_hpp
#define Vec3dVector_hpp

#include <stddef.h>
#include <stdbool.h>

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

typedef char* Vec3dVector_ptr;

#ifdef __cplusplus
extern "C" {
#endif

    DLLEXPORT Vec3dVector_ptr Vec3dVector_new(size_t s);
    DLLEXPORT double* Vec3dVector_get(Vec3dVector_ptr self, size_t idx);
    DLLEXPORT size_t Vec3dVector_size(Vec3dVector_ptr self);
    DLLEXPORT void Vec3dVector_delete(Vec3dVector_ptr self);

#ifdef __cplusplus
}
#endif
#endif /* Vec3dVector_hpp */
