//
//  Vec3dVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Vec3dVector_hpp
#define Vec3dVector_hpp

#include <vector>
#include <cstddef>
#include <GEL/CGLA/Vec3d.h>

namespace PyGEL {
    using namespace CGLA;
    using Vec3dVector = std::vector<Vec3d>;
    using Vec3dVector_ptr = Vec3dVector*; // C-style alias
    
    Vec3dVector_ptr Vec3dVector_new(size_t s);
    double* Vec3dVector_get(Vec3dVector_ptr self, size_t idx);
    size_t Vec3dVector_size(Vec3dVector_ptr self);
    void Vec3dVector_delete(Vec3dVector_ptr self);
}

#endif /* Vec3dVector_hpp */
