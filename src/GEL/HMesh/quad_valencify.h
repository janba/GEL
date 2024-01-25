//
//  comb_quad.hpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 05/03/2023.
//  Copyright © 2023 J. Andreas Bærentzen. All rights reserved.
//

#ifndef comb_quad_hpp
#define comb_quad_hpp

#include <stdio.h>

#endif /* comb_quad_hpp */

#include <GEL/HMesh/HMesh.h>

/** This function turns the input Manifold m into a mesh where all vertices have valence 4.
 This is typically not in itself very interesting, but the dual  mesh is then a quad mesh, and that
 is useful for inverse skeletonization. */
void quad_valencify(HMesh::Manifold& m);