//
//  PyGEL.c
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 02/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "I3DTree.h"
#include "Vec3dVector.h"


#include <GEL/CGLA/Vec3d.h>
#include <GEL/Geometry/KDTree.h>

using namespace PyGEL;
using namespace CGLA;
using namespace Geometry;

using namespace CGLA;
using namespace std;

template class Geometry::KDTree<CGLA::Vec3d, size_t>;

namespace PyGEL {
    using Vec3dVector = vector<Vec3d>;


    I3DTree_ptr I3DTree_new() {
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


    std::pair<std::vector<double>, size_t> I3DTree_closest_point(I3DTree* tree, double x, double y, double z, double r) {
        CGLA::Vec3d p(x, y, z);
        CGLA::Vec3d key;
        size_t val = 0;
        if (tree->closest_point(p, r, key, val)) {
            return { {key[0], key[1], key[2]}, val };
        }
        return { {}, 0 };
    }


    std::pair<Vec3dVector, std::vector<size_t>> I3DTree_in_sphere(I3DTree* tree, double x, double y, double z, double r) {
        Vec3dVector keys;
        std::vector<size_t> vals;
        tree->in_sphere(Vec3d(x,y,z), r, keys, vals);
        return {keys, vals};
    }


    std::pair<Vec3dVector, std::vector<size_t>> I3DTree_m_closest_points(I3DTree* tree, double x, double y, double z, double r, int m) {
        CGLA::Vec3d p(x,y,z);
        auto records = tree->m_closest(m,p,r);
        Vec3dVector keys(records.size());
        std::vector<size_t> vals(records.size());
        for(size_t i=0; i<records.size(); ++i) {
            keys[i] = records[i].k;
            vals[i] = records[i].v;
        }
        return {keys, vals};
    }
}