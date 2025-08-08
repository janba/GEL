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

#include <GEL/CGLA/Vec3d.h>
#include <GEL/Geometry/KDTree.h>

using namespace PyGEL;
using namespace CGLA;
using namespace Geometry;

using namespace CGLA;
using namespace std;

template class Geometry::KDTree<CGLA::Vec3d, size_t>;
using Vec3dVector = vector<Vec3d>;
using IntVector = vector<size_t>;

I3DTree_ptr I3DTree_new() {
    I3DTree* ptr = new I3DTree();
    return ptr;
}

void I3DTree_delete(I3DTree_ptr self) {
    delete reinterpret_cast<I3DTree*>(self);
}

void I3DTree_insert(I3DTree_ptr tree, double x, double y, double z, size_t v) {
    tree->insert(CGLA::Vec3d(x,y,z), v);
}

void I3DTree_build(I3DTree_ptr tree) {
    tree->build();
}

size_t I3DTree_closest_point(I3DTree_ptr tree, double x, double y, double z, double r,
                                  double* key, size_t* val) {
    CGLA::Vec3d p(x,y,z);
    if(tree->closest_point(p, r,
                                                       *reinterpret_cast<Vec3d*>(key),
                                                       *val))
        return 1;
    return 0;
}

size_t I3DTree_in_sphere(I3DTree_ptr tree, double x, double y, double z, double r,
                              Vec3dVector_ptr keys, IntVector_ptr vals) {
    return tree->in_sphere(Vec3d(x,y,z), r,
                                                       *reinterpret_cast<Vec3dVector*>(keys),
                                                       *reinterpret_cast<IntVector*>(vals));
}

size_t I3DTree_m_closest_points(I3DTree_ptr tree, double x, double y, double z, double r, int m,
                                Vec3dVector_ptr _keys, IntVector_ptr _vals){
    CGLA::Vec3d p(x,y,z);
    auto records = tree->m_closest(m,p,r);
    
    Vec3dVector& keys = *reinterpret_cast<Vec3dVector*>(_keys);
    IntVector& vals = *reinterpret_cast<IntVector*>(_vals);

    auto N = records.size();
    keys.resize(N);
    vals.resize(N);
    for(int i=0;i<records.size();++i) {
        keys[i] = records[i].k;
        vals[i] = records[i].v;
    }
    return N;
}
