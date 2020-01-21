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

#include <stdbool.h>
#include "IntVector.h"
#include "Vec3dVector.h"

typedef char* Manifold_ptr;

#ifdef __cplusplus
extern "C" {
#endif

//    Manifold class methods
    DLLEXPORT Manifold_ptr Manifold_new();

    DLLEXPORT Manifold_ptr Manifold_from_triangles(int NV, int NF, double* vertices, int* faces);

    DLLEXPORT Manifold_ptr Manifold_from_points(int N, double* pts, double* X_axis, double* Y_axis);

    DLLEXPORT Manifold_ptr Manifold_copy(Manifold_ptr self);

    DLLEXPORT void Manifold_delete(Manifold_ptr self);
    DLLEXPORT size_t Manifold_positions(Manifold_ptr self, double** pos);
    
    DLLEXPORT size_t Manifold_no_allocated_vertices(Manifold_ptr self);
    DLLEXPORT size_t Manifold_no_allocated_faces(Manifold_ptr self);
    DLLEXPORT size_t Manifold_no_allocated_halfedges(Manifold_ptr self);
    DLLEXPORT size_t Manifold_vertices(Manifold_ptr self, IntVector_ptr verts);
    DLLEXPORT size_t Manifold_faces(Manifold_ptr self, IntVector_ptr faces);
    DLLEXPORT size_t Manifold_halfedges(Manifold_ptr self, IntVector_ptr hedges);
    DLLEXPORT size_t Manifold_circulate_vertex(Manifold_ptr self, size_t v, char mode, IntVector_ptr nverts);
    DLLEXPORT size_t Manifold_circulate_face(Manifold_ptr self, size_t f, char mode, IntVector_ptr nverts);
    
    DLLEXPORT void Manifold_add_face(Manifold_ptr self, size_t no_verts, double* pos);

    DLLEXPORT bool Manifold_remove_face(Manifold_ptr self,size_t fid);
    DLLEXPORT bool Manifold_remove_edge(Manifold_ptr self,size_t hid);
    DLLEXPORT bool Manifold_remove_vertex(Manifold_ptr self,size_t vid);

    DLLEXPORT bool Manifold_vertex_in_use(Manifold_ptr self,size_t id);
    DLLEXPORT bool Manifold_face_in_use(Manifold_ptr self,size_t id);
    DLLEXPORT bool Manifold_halfedge_in_use(Manifold_ptr self,size_t id);
    
    DLLEXPORT bool Manifold_flip_edge(Manifold_ptr self,size_t h);
    DLLEXPORT bool Manifold_collapse_edge(Manifold_ptr self,size_t h, bool avg_vertices);

    DLLEXPORT size_t Manifold_split_face_by_edge(Manifold_ptr self,size_t f, size_t v0, size_t v1);
    DLLEXPORT size_t Manifold_split_face_by_vertex(Manifold_ptr self,size_t f);
    DLLEXPORT size_t Manifold_split_edge(Manifold_ptr self,size_t h);
    
    DLLEXPORT bool Manifold_stitch_boundary_edges(Manifold_ptr self,size_t h0, size_t h1);
    DLLEXPORT bool Manifold_merge_faces(Manifold_ptr self,size_t f, size_t h);

    DLLEXPORT size_t Manifold_close_hole(Manifold_ptr self,size_t h);
    DLLEXPORT void Manifold_cleanup(Manifold_ptr self);
    
//  New functions that rely on walker
    
    DLLEXPORT size_t Walker_next_halfedge(Manifold_ptr m_ptr, size_t h);
    DLLEXPORT size_t Walker_prev_halfedge(Manifold_ptr m_ptr, size_t h);
    DLLEXPORT size_t Walker_opposite_halfedge(Manifold_ptr m_ptr, size_t h);
    DLLEXPORT size_t Walker_incident_face(Manifold_ptr m_ptr, size_t h);
    DLLEXPORT size_t Walker_incident_vertex(Manifold_ptr m_ptr, size_t h);
    
// External functions
    
    DLLEXPORT bool is_halfedge_at_boundary(const Manifold_ptr m_ptr, size_t _h);
    DLLEXPORT bool is_vertex_at_boundary(const Manifold_ptr m_ptr, size_t _v);
    
    DLLEXPORT double length(const Manifold_ptr m_ptr, size_t _h);
    DLLEXPORT bool boundary_edge(const Manifold_ptr m_ptr, size_t _v, size_t _h);
    DLLEXPORT size_t valency(const Manifold_ptr m_ptr, size_t _v);
    DLLEXPORT void vertex_normal(const Manifold_ptr m_ptr, size_t _v, double*);
    DLLEXPORT bool connected(const Manifold_ptr m_ptr, size_t _v0, size_t _v1);
    
    DLLEXPORT size_t no_edges(const Manifold_ptr m_ptr, size_t _f);
    DLLEXPORT void face_normal(const Manifold_ptr m_ptr, size_t _f, double*);
    DLLEXPORT double area(const Manifold_ptr m_ptr, size_t _f);
    DLLEXPORT double perimeter(const Manifold_ptr m_ptr, size_t _f);
    DLLEXPORT void centre(const Manifold_ptr m_ptr, size_t _f, double*);
    
    DLLEXPORT extern size_t InvalidIndex;

#ifdef __cplusplus
}
#endif

#endif /* Manifold_hpp */
