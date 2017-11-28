//
//  hmesh_functions.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "hmesh_functions.h"
#include <string>

using namespace std;

void obj_load(char* fn, HMesh::Manifold* m_ptr) {
    obj_load(string(fn), *m_ptr, true);
}
