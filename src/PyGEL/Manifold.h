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
#include "IntVector.h"
#include "Vec3dVector.h"

namespace PyGEL {
    namespace py = pybind11;
    using namespace HMesh;
    using Vec = HMesh::Manifold::Vec;
    using Scalar = HMesh::Manifold::Vec::ScalarType;
    using Manifold_ptr = Manifold*;
    using IntVector_ptr = IntVector*; // C-style alias
    
    // Manifold class methods
    Manifold_ptr Manifold_new();
    Manifold_ptr Manifold_from_triangles(const std::vector<double>& vertices, const std::vector<int>& faces);
    Manifold_ptr Manifold_from_points(int N, const std::vector<double>& pts, const Vec& X_axis, const Vec& Y_axis);
    Manifold_ptr Manifold_copy(Manifold_ptr self);
    void Manifold_merge(Manifold_ptr self, Manifold_ptr other);
    void Manifold_delete(Manifold_ptr self);
    py::array_t<Scalar> Manifold_positions(Manifold_ptr self);
    
    size_t Manifold_no_allocated_vertices(Manifold_ptr self);
    size_t Manifold_no_allocated_faces(Manifold_ptr self);
    size_t Manifold_no_allocated_halfedges(Manifold_ptr self);
    size_t Manifold_vertices(Manifold_ptr self, IntVector& verts);
    size_t Manifold_faces(Manifold_ptr self, IntVector& faces);
    size_t Manifold_halfedges(Manifold_ptr self, IntVector& hedges);
    size_t Manifold_circulate_vertex(Manifold_ptr self, size_t v, char mode, IntVector& nverts);
    size_t Manifold_circulate_face(Manifold_ptr self, size_t f, char mode, IntVector& nverts);
    
    size_t Manifold_add_face(Manifold_ptr self, const std::vector<double>& pos);
    bool Manifold_remove_face(Manifold_ptr self, size_t fid);
    bool Manifold_remove_edge(Manifold_ptr self, size_t hid);
    bool Manifold_remove_vertex(Manifold_ptr self, size_t vid);
    bool Manifold_vertex_in_use(Manifold_ptr self, size_t id);
    bool Manifold_face_in_use(Manifold_ptr self, size_t id);
    bool Manifold_halfedge_in_use(Manifold_ptr self, size_t id);
    
    bool Manifold_flip_edge(Manifold_ptr self, size_t h);
    bool Manifold_collapse_edge(Manifold_ptr self, size_t h, bool avg_vertices);
    size_t Manifold_split_face_by_edge(Manifold_ptr self, size_t f, size_t v0, size_t v1);
    size_t Manifold_split_face_by_vertex(Manifold_ptr self, size_t f);
    size_t Manifold_split_edge(Manifold_ptr self, size_t h);
    
    bool Manifold_stitch_boundary_edges(Manifold_ptr self, size_t h0, size_t h1);
    bool Manifold_merge_faces(Manifold_ptr self, size_t f, size_t h);
    size_t Manifold_close_hole(Manifold_ptr self, size_t h);
    void Manifold_cleanup(Manifold_ptr self);
    
    // Walker functions
    size_t Walker_next_halfedge(Manifold_ptr m_ptr, size_t h);
    size_t Walker_prev_halfedge(Manifold_ptr m_ptr, size_t h);
    size_t Walker_opposite_halfedge(Manifold_ptr m_ptr, size_t h);
    size_t Walker_incident_face(Manifold_ptr m_ptr, size_t h);
    size_t Walker_incident_vertex(Manifold_ptr m_ptr, size_t h);
    
    // External functions
    bool is_halfedge_at_boundary(const Manifold_ptr m_ptr, size_t _h);
    bool is_vertex_at_boundary(const Manifold_ptr m_ptr, size_t _v);
    double length(const Manifold_ptr m_ptr, size_t _h);
    bool boundary_edge(const Manifold_ptr m_ptr, size_t _v, size_t _h);
    size_t valency(const Manifold_ptr m_ptr, size_t _v);
    Vec vertex_normal(const Manifold_ptr _m_ptr, size_t _v);
    bool connected(const Manifold_ptr m_ptr, size_t _v0, size_t _v1);
    
    size_t no_edges(const Manifold_ptr m_ptr, size_t _f);
    Vec face_normal(const Manifold_ptr m_ptr, size_t _f);
    double area(const Manifold_ptr m_ptr, size_t _f);
    double one_ring_area(const Manifold_ptr m_ptr, size_t _v);
    double mixed_area(const Manifold_ptr m_ptr, size_t _v);
    double gaussian_curvature(const Manifold_ptr m_ptr, size_t _v);
    double mean_curvature(const Manifold_ptr m_ptr, size_t _v);
    std::vector<double> principal_curvatures(const Manifold_ptr m_ptr, size_t _v);
    double perimeter(const Manifold_ptr m_ptr, size_t _f);
    Vec centre(const Manifold_ptr m_ptr, size_t _f);
    double total_area(const Manifold_ptr);
    double volume(const Manifold_ptr);
    
    extern size_t InvalidIndex;
}

#endif /* Manifold_hpp */
