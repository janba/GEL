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

void quad_valencify_no_crossings(HMesh::Manifold& m,
                                 HMesh::VertexAttributeVector<CGLA::Vec3f>&,
                                 HMesh::HalfEdgeAttributeVector<CGLA::Vec3f>&);

void quad_valencify_cc(HMesh::Manifold& m_orig,
                       HMesh::VertexAttributeVector<CGLA::Vec3f>& vcol,
                       HMesh::HalfEdgeAttributeVector<CGLA::Vec3f>& hcol,
                       HMesh::FaceAttributeVector<CGLA::Vec3f>& fcol);
