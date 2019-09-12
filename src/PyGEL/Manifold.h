//
//  Manifold.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 11/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Manifold_hpp
#define Manifold_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif
    
#include <GEL/HMesh/HMesh.h>
#include "IntVector.h"
#include "Vec3dVector.h"

extern "C" {
    
//    Manifold class methods
    DLLEXPORT HMesh::Manifold* Manifold_new();

    DLLEXPORT HMesh::Manifold* Manifold_from_triangles(int NV, int NF, double* vertices, int* faces);

    DLLEXPORT HMesh::Manifold* Manifold_from_points(int N, double* pts, double* X_axis, double* Y_axis);

    DLLEXPORT HMesh::Manifold* Manifold_copy(HMesh::Manifold* self);

    DLLEXPORT void Manifold_delete(HMesh::Manifold* self);
    DLLEXPORT size_t Manifold_positions(HMesh::Manifold* self, double** pos);
    
    DLLEXPORT size_t Manifold_no_allocated_vertices(HMesh::Manifold* self);
    DLLEXPORT size_t Manifold_no_allocated_faces(HMesh::Manifold* self);
    DLLEXPORT size_t Manifold_no_allocated_halfedges(HMesh::Manifold* self);
    DLLEXPORT size_t Manifold_vertices(HMesh::Manifold* self, IntVector* verts);
    DLLEXPORT size_t Manifold_faces(HMesh::Manifold* self, IntVector* faces);
    DLLEXPORT size_t Manifold_halfedges(HMesh::Manifold* self, IntVector* hedges);
    DLLEXPORT size_t Manifold_circulate_vertex(HMesh::Manifold* self, size_t v, char mode, IntVector* nverts);
    DLLEXPORT size_t Manifold_circulate_face(HMesh::Manifold* self, size_t f, char mode, IntVector* nverts);
    
    DLLEXPORT void Manifold_add_face(HMesh::Manifold* self, size_t no_verts, double* pos);

    DLLEXPORT bool Manifold_remove_face(HMesh::Manifold* self,size_t fid);
    DLLEXPORT bool Manifold_remove_edge(HMesh::Manifold* self,size_t hid);
    DLLEXPORT bool Manifold_remove_vertex(HMesh::Manifold* self,size_t vid);

    DLLEXPORT bool Manifold_vertex_in_use(HMesh::Manifold* self,size_t id);
    DLLEXPORT bool Manifold_face_in_use(HMesh::Manifold* self,size_t id);
    DLLEXPORT bool Manifold_halfedge_in_use(HMesh::Manifold* self,size_t id);
    
    DLLEXPORT bool Manifold_flip_edge(HMesh::Manifold* self,size_t h);
    DLLEXPORT bool Manifold_collapse_edge(HMesh::Manifold* self,size_t h, bool avg_vertices = false);

    DLLEXPORT size_t Manifold_split_face_by_edge(HMesh::Manifold* self,size_t f, size_t v0, size_t v1);
    DLLEXPORT size_t Manifold_split_face_by_vertex(HMesh::Manifold* self,size_t f);
    DLLEXPORT size_t Manifold_split_edge(HMesh::Manifold* self,size_t h);
    
    DLLEXPORT bool Manifold_stitch_boundary_edges(HMesh::Manifold* self,size_t h0, size_t h1);
    DLLEXPORT bool Manifold_merge_faces(HMesh::Manifold* self,size_t f, size_t h);

    DLLEXPORT size_t Manifold_close_hole(HMesh::Manifold* self,size_t h);
    DLLEXPORT void Manifold_cleanup(HMesh::Manifold* self);
    
//  New functions that rely on walker
    
    DLLEXPORT size_t Walker_next_halfedge(HMesh::Manifold* m_ptr, size_t h);
    DLLEXPORT size_t Walker_prev_halfedge(HMesh::Manifold* m_ptr, size_t h);
    DLLEXPORT size_t Walker_opposite_halfedge(HMesh::Manifold* m_ptr, size_t h);
    DLLEXPORT size_t Walker_incident_face(HMesh::Manifold* m_ptr, size_t h);
    DLLEXPORT size_t Walker_incident_vertex(HMesh::Manifold* m_ptr, size_t h);
    
// External functions
    
    DLLEXPORT bool is_halfedge_at_boundary(const HMesh::Manifold* m_ptr, size_t _h);
    DLLEXPORT bool is_vertex_at_boundary(const HMesh::Manifold* m_ptr, size_t _v);
    
    DLLEXPORT double length(const HMesh::Manifold* m_ptr, size_t _h);
    DLLEXPORT bool boundary_edge(const HMesh::Manifold* m_ptr, size_t _v, size_t _h);
    DLLEXPORT size_t valency(const HMesh::Manifold* m_ptr, size_t _v);
    DLLEXPORT void vertex_normal(const HMesh::Manifold* m_ptr, size_t _v, CGLA::Vec3d*);
    DLLEXPORT bool connected(const HMesh::Manifold* m_ptr, size_t _v0, size_t _v1);
    
    DLLEXPORT size_t no_edges(const HMesh::Manifold* m_ptr, size_t _f);
    DLLEXPORT void face_normal(const HMesh::Manifold* m_ptr, size_t _f, CGLA::Vec3d*);
    DLLEXPORT double area(const HMesh::Manifold* m_ptr, size_t _f);
    DLLEXPORT double perimeter(const HMesh::Manifold* m_ptr, size_t _f);
    DLLEXPORT void centre(const HMesh::Manifold* m_ptr, size_t _f, CGLA::Vec3d*);
    
    DLLEXPORT extern size_t InvalidIndex;
}

#endif /* Manifold_hpp */
