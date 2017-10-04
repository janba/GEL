//
//  Vec3dVector.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "Vec3dVector.h"
#include <iostream>
using namespace CGLA;
using namespace std;

Vec3dVector* Vec3dVector_new(size_t s) {
    return new Vec3dVector(s);
}

double* Vec3dVector_get(Vec3dVector* self, int idx) {
    return (*self)[idx].get();
}

void Vec3dVector_delete(Vec3dVector* self) {
    delete self;
}
