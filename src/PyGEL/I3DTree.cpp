//
//  PyGEL.c
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 02/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "I3DTree.h"
#include "Vec3dVector.h"
#include "IntVector.h"
#include <iostream>
using namespace CGLA;
using namespace std;

I3DTree* I3DTree_new() {
    I3DTree* ptr = new I3DTree();
    return ptr;
}
void I3DTree_delete(I3DTree* self) {
    delete self;
}

void I3DTree_insert(I3DTree* tree, double x, double y, double z, size_t v) {
    tree->insert(CGLA::Vec3d(x,y,z), v);
}

void I3DTree_build(I3DTree* tree) {
    tree->build();
}

size_t I3DTree_closest_point(I3DTree* tree, double x, double y, double z, double r,
                             CGLA::Vec3d* key, size_t* val) {
    CGLA::Vec3d p(x,y,z);
    if(tree->closest_point(p, r, *key, *val))
        return 1;
    return 0;
}

size_t I3DTree_in_sphere(I3DTree* tree, double x, double y, double z, double r,
                         Vec3dVector* keys, IntVector* vals) {
    return tree->in_sphere(Vec3d(x,y,z), r, *keys, *vals);
}

size_t I3DTree_m_closest_points(I3DTree* tree, double x, double y, double z, double r, int m,
                                CGLA::Vec3d* key, size_t* val){
    CGLA::Vec3d p(x,y,z);
    auto records = tree->m_closest(m,p,r);
        return 1;
    return 0;
}
