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

    Vec3dVector* Vec3dVector_new(size_t s);
    double* Vec3dVector_get(Vec3dVector* self, size_t idx);
    size_t Vec3dVector_size(Vec3dVector* self);
    void Vec3dVector_delete(Vec3dVector* self);
}

#endif /* Vec3dVector_hpp */
