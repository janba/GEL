//
//  Manifold.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 11/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Manifold_hpp
#define Manifold_hpp

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <GEL/HMesh/Manifold.h>

#include "Vec3dVector.h"

namespace PyGEL {
    namespace py = pybind11;
    using Vec = HMesh::Manifold::Vec;
    using Scalar = HMesh::Manifold::Vec::ScalarType;


    
    // Manifold class methods
    // HMesh::Manifold* Manifold_new();
    HMesh::Manifold* Manifold_from_triangles(const std::vector<double>& vertices, const std::vector<int>& faces);
    HMesh::Manifold* Manifold_from_points(int N, const std::vector<double>& pts, const Vec& X_axis, const Vec& Y_axis);
    // HMesh::Manifold* Manifold_copy(HMesh::Manifold* self);
    void Manifold_merge(HMesh::Manifold* self, HMesh::Manifold* other);
    void Manifold_delete(HMesh::Manifold* self);
    py::array_t<Scalar> Manifold_positions(HMesh::Manifold* self);

    size_t Manifold_no_allocated_vertices(HMesh::Manifold* self);
    size_t Manifold_no_allocated_faces(HMesh::Manifold* self);
    size_t Manifold_no_allocated_halfedges(HMesh::Manifold* self);
    std::vector<size_t> Manifold_vertices(HMesh::Manifold* self);
    std::vector<size_t> Manifold_faces(HMesh::Manifold* self);
    std::vector<size_t> Manifold_halfedges(HMesh::Manifold* self);
    std::vector<size_t> Manifold_circulate_vertex(HMesh::Manifold* self, size_t v, char mode);
    std::vector<size_t> Manifold_circulate_face(HMesh::Manifold* self, size_t f, char mode);

    size_t Manifold_add_face(HMesh::Manifold* self, const std::vector<double>& pos);
    bool Manifold_remove_face(HMesh::Manifold* self, size_t fid);
    bool Manifold_remove_edge(HMesh::Manifold* self, size_t hid);
    bool Manifold_remove_vertex(HMesh::Manifold* self, size_t vid);
    bool Manifold_vertex_in_use(HMesh::Manifold* self, size_t id);
    bool Manifold_face_in_use(HMesh::Manifold* self, size_t id);
    bool Manifold_halfedge_in_use(HMesh::Manifold* self, size_t id);

    bool Manifold_flip_edge(HMesh::Manifold* self, size_t h);
    bool Manifold_collapse_edge(HMesh::Manifold* self, size_t h, bool avg_vertices);
    size_t Manifold_split_face_by_edge(HMesh::Manifold* self, size_t f, size_t v0, size_t v1);
    size_t Manifold_split_face_by_vertex(HMesh::Manifold* self, size_t f);
    size_t Manifold_split_edge(HMesh::Manifold* self, size_t h);

    bool Manifold_stitch_boundary_edges(HMesh::Manifold* self, size_t h0, size_t h1);
    bool Manifold_merge_faces(HMesh::Manifold* self, size_t f, size_t h);
    size_t Manifold_close_hole(HMesh::Manifold* self, size_t h);
    void Manifold_cleanup(HMesh::Manifold* self);
    
    // Walker functions
    size_t Walker_next_halfedge(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_prev_halfedge(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_opposite_halfedge(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_incident_face(HMesh::Manifold* m_ptr, size_t h);
    size_t Walker_incident_vertex(HMesh::Manifold* m_ptr, size_t h);
    
    // External functions
    bool is_halfedge_at_boundary(const HMesh::Manifold* m_ptr, size_t _h);
    bool is_vertex_at_boundary(const HMesh::Manifold* m_ptr, size_t _v);
    double length(const HMesh::Manifold* m_ptr, size_t _h);
    bool boundary_edge(const HMesh::Manifold* m_ptr, size_t _v, size_t _h);
    size_t valency(const HMesh::Manifold* m_ptr, size_t _v);
    Vec vertex_normal(const HMesh::Manifold* _m_ptr, size_t _v);
    bool connected(const HMesh::Manifold* m_ptr, size_t _v0, size_t _v1);
    
    size_t no_edges(const HMesh::Manifold* m_ptr, size_t _f);
    Vec face_normal(const HMesh::Manifold* m_ptr, size_t _f);
    double area(const HMesh::Manifold* m_ptr, size_t _f);
    double one_ring_area(const HMesh::Manifold* m_ptr, size_t _v);
    double mixed_area(const HMesh::Manifold* m_ptr, size_t _v);
    double gaussian_curvature(const HMesh::Manifold* m_ptr, size_t _v);
    double mean_curvature(const HMesh::Manifold* m_ptr, size_t _v);
    std::vector<double> principal_curvatures(const HMesh::Manifold* m_ptr, size_t _v);
    double perimeter(const HMesh::Manifold* m_ptr, size_t _f);
    Vec centre(const HMesh::Manifold* m_ptr, size_t _f);
    double total_area(const HMesh::Manifold*);
    double volume(const HMesh::Manifold*);
    
    extern size_t InvalidIndex;
}

#endif /* Manifold_hpp */
