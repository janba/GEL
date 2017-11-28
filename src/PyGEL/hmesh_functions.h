//
//  hmesh_functions.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef hmesh_functions_hpp
#define hmesh_functions_hpp

#include <GEL/HMesh/HMesh.h>

extern "C" {
    void obj_load(char*, HMesh::Manifold*);
}

#endif /* hmesh_functions_hpp */
