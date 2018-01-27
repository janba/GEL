//
//  Vec3dVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Vec3dVector_hpp
#define Vec3dVector_hpp

#ifdef __APPLE__
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <GEL/CGLA/Vec3d.h>
#include <vector>

using Vec3dVector = std::vector<CGLA::Vec3d>;

extern "C" {
    DLLEXPORT Vec3dVector* Vec3dVector_new(size_t s);
    DLLEXPORT double* Vec3dVector_get(Vec3dVector* self, size_t idx);
    DLLEXPORT size_t Vec3dVector_size(Vec3dVector* self);
    DLLEXPORT void Vec3dVector_delete(Vec3dVector* self);
}
#endif /* Vec3dVector_hpp */
