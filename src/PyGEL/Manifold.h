//
//  Manifold.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 11/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Manifold_hpp
#define Manifold_hpp

#include <GEL/HMesh/HMesh.h>
#include "IntVector.h"
#include "Vec3dVector.h"

extern "C" {
    
//    Manifold class methods
    HMesh::Manifold* Manifold_new();
    
    HMesh::Manifold* Manifold_copy(HMesh::Manifold* self);

    void Manifold_delete(HMesh::Manifold* self);
    size_t Manifold_positions(HMesh::Manifold* self, double** pos);
    
    size_t Manifold_no_allocated_vertices(HMesh::Manifold* self);
    size_t Manifold_no_allocated_faces(HMesh::Manifold* self);
    size_t Manifold_no_allocated_halfedges(HMesh::Manifold* self);
    size_t Manifold_vertices(HMesh::Manifold* self, IntVector* verts);
    size_t Manifold_faces(HMesh::Manifold* self, IntVector* faces);
    size_t Manifold_halfedges(HMesh::Manifold* self, IntVector* hedges);
    size_t Manifold_circulate_vertex(HMesh::Manifold* self, size_t v, char mode, IntVector* nverts);
    size_t Manifold_circulate_face(HMesh::Manifold* self, size_t f, char mode, IntVector* nverts);
    
    void Manifold_add_face(HMesh::Manifold* self, size_t no_verts, double* pos);

    bool Manifold_remove_face(HMesh::Manifold* self,size_t fid);
    bool Manifold_remove_edge(HMesh::Manifold* self,size_t hid);
    bool Manifold_remove_vertex(HMesh::Manifold* self,size_t vid);

    bool Manifold_vertex_in_use(HMesh::Manifold* self,size_t id);
    bool Manifold_face_in_use(HMesh::Manifold* self,size_t id);
    bool Manifold_halfedge_in_use(HMesh::Manifold* self,size_t id);
    
    bool Manifold_flip_edge(HMesh::Manifold* self,size_t h);
    bool Manifold_collapse_edge(HMesh::Manifold* self,size_t h, bool avg_vertices = false);

    size_t Manifold_split_face_by_edge(HMesh::Manifold* self,size_t f, size_t v0, size_t v1);
    size_t Manifold_split_face_by_vertex(HMesh::Manifold* self,size_t f);
    size_t Manifold_split_edge(HMesh::Manifold* self,size_t h);
    
    bool Manifold_stitch_boundary_edges(HMesh::Manifold* self,size_t h0, size_t h1);
    bool Manifold_merge_faces(HMesh::Manifold* self,size_t f, size_t h);

    size_t Manifold_close_hole(HMesh::Manifold* self,size_t h);
    void Manifold_cleanup(HMesh::Manifold* self);
    
//  New functions that rely on walker
    
    size_t Walker_next_halfedge(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_opposite_halfedge(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_incident_face(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_incident_vertex(HMesh::Manifold* m_ptr, size_t h);
    
// External functions
    
    bool is_halfedge_at_boundary(const HMesh::Manifold* m_ptr, size_t _h);
    bool is_vertex_at_boundary(const HMesh::Manifold* m_ptr, size_t _v);
    
    double length(const HMesh::Manifold* m_ptr, size_t _h);
    bool boundary_edge(const HMesh::Manifold* m_ptr, size_t _v, size_t _h);
    size_t valency(const HMesh::Manifold* m_ptr, size_t _v);
    void vertex_normal(const HMesh::Manifold* m_ptr, size_t _v, CGLA::Vec3d*);
    bool connected(const HMesh::Manifold* m_ptr, size_t _v0, size_t _v1);
    
    size_t no_edges(const HMesh::Manifold* m_ptr, size_t _f);
    void face_normal(const HMesh::Manifold* m_ptr, size_t _f, CGLA::Vec3d*);
    double area(const HMesh::Manifold* m_ptr, size_t _f);
    double perimeter(const HMesh::Manifold* m_ptr, size_t _f);
    void centre(const HMesh::Manifold* m_ptr, size_t _f, CGLA::Vec3d*);
        
}

#endif /* Manifold_hpp */
