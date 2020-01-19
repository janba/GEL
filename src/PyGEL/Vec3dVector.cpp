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

using Vec3dVector = vector<Vec3d>;

Vec3dVector_ptr Vec3dVector_new(size_t s) {
    return reinterpret_cast<Vec3dVector_ptr>(new Vec3dVector(s));
}

size_t Vec3dVector_size(Vec3dVector_ptr self) {
    return reinterpret_cast<Vec3dVector*>(self)->size();
}


double* Vec3dVector_get(Vec3dVector_ptr self, size_t idx) {
    return (*reinterpret_cast<Vec3dVector*>(self))[idx].get();
}

void Vec3dVector_delete(Vec3dVector_ptr self) {
    delete reinterpret_cast<Vec3dVector*>(self);
}
