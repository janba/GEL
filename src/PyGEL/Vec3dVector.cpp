//
//  Vec3dVector.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "Vec3dVector.h"

#include <GEL/CGLA/CGLA.h>
#include <vector>

using namespace CGLA;
using namespace std;

namespace PyGEL {

Vec3dVector* Vec3dVector_new(size_t s) {
    return new Vec3dVector(s);
}

size_t Vec3dVector_size(Vec3dVector* self) {
    return self->size();
}

double* Vec3dVector_get(Vec3dVector* self, size_t idx) {
    return (*self)[idx].get();
}

void Vec3dVector_delete(Vec3dVector* self) {
    delete self;
}

} // namespace PyGEL
