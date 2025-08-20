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
    using namespace HMesh;
    using Vec = HMesh::Manifold::Vec;
    using Scalar = HMesh::Manifold::Vec::ScalarType;


    
    // Manifold class methods
    Manifold* Manifold_new();
    Manifold* Manifold_from_triangles(const std::vector<double>& vertices, const std::vector<int>& faces);
    Manifold* Manifold_from_points(int N, const std::vector<double>& pts, const Vec& X_axis, const Vec& Y_axis);
    Manifold* Manifold_copy(Manifold* self);
    void Manifold_merge(Manifold* self, Manifold* other);
    void Manifold_delete(Manifold* self);
    py::array_t<Scalar> Manifold_positions(Manifold* self);

    size_t Manifold_no_allocated_vertices(Manifold* self);
    size_t Manifold_no_allocated_faces(Manifold* self);
    size_t Manifold_no_allocated_halfedges(Manifold* self);
    std::vector<size_t> Manifold_vertices(Manifold* self);
    std::vector<size_t> Manifold_faces(Manifold* self);
    std::vector<size_t> Manifold_halfedges(Manifold* self);
    std::vector<size_t> Manifold_circulate_vertex(Manifold* self, size_t v, char mode);
    std::vector<size_t> Manifold_circulate_face(Manifold* self, size_t f, char mode);

    size_t Manifold_add_face(Manifold* self, const std::vector<double>& pos);
    bool Manifold_remove_face(Manifold* self, size_t fid);
    bool Manifold_remove_edge(Manifold* self, size_t hid);
    bool Manifold_remove_vertex(Manifold* self, size_t vid);
    bool Manifold_vertex_in_use(Manifold* self, size_t id);
    bool Manifold_face_in_use(Manifold* self, size_t id);
    bool Manifold_halfedge_in_use(Manifold* self, size_t id);

    bool Manifold_flip_edge(Manifold* self, size_t h);
    bool Manifold_collapse_edge(Manifold* self, size_t h, bool avg_vertices);
    size_t Manifold_split_face_by_edge(Manifold* self, size_t f, size_t v0, size_t v1);
    size_t Manifold_split_face_by_vertex(Manifold* self, size_t f);
    size_t Manifold_split_edge(Manifold* self, size_t h);

    bool Manifold_stitch_boundary_edges(Manifold* self, size_t h0, size_t h1);
    bool Manifold_merge_faces(Manifold* self, size_t f, size_t h);
    size_t Manifold_close_hole(Manifold* self, size_t h);
    void Manifold_cleanup(Manifold* self);
    
    // Walker functions
    size_t Walker_next_halfedge(Manifold* m_ptr, size_t h);
    size_t Walker_prev_halfedge(Manifold* m_ptr, size_t h);
    size_t Walker_opposite_halfedge(Manifold* m_ptr, size_t h);
    size_t Walker_incident_face(Manifold* m_ptr, size_t h);
    size_t Walker_incident_vertex(Manifold* m_ptr, size_t h);
    
    // External functions
    bool is_halfedge_at_boundary(const Manifold* m_ptr, size_t _h);
    bool is_vertex_at_boundary(const Manifold* m_ptr, size_t _v);
    double length(const Manifold* m_ptr, size_t _h);
    bool boundary_edge(const Manifold* m_ptr, size_t _v, size_t _h);
    size_t valency(const Manifold* m_ptr, size_t _v);
    Vec vertex_normal(const Manifold* _m_ptr, size_t _v);
    bool connected(const Manifold* m_ptr, size_t _v0, size_t _v1);
    
    size_t no_edges(const Manifold* m_ptr, size_t _f);
    Vec face_normal(const Manifold* m_ptr, size_t _f);
    double area(const Manifold* m_ptr, size_t _f);
    double one_ring_area(const Manifold* m_ptr, size_t _v);
    double mixed_area(const Manifold* m_ptr, size_t _v);
    double gaussian_curvature(const Manifold* m_ptr, size_t _v);
    double mean_curvature(const Manifold* m_ptr, size_t _v);
    std::vector<double> principal_curvatures(const Manifold* m_ptr, size_t _v);
    double perimeter(const Manifold* m_ptr, size_t _f);
    Vec centre(const Manifold* m_ptr, size_t _f);
    double total_area(const Manifold*);
    double volume(const Manifold*);
    
    extern size_t InvalidIndex;
}

#endif /* Manifold_hpp */
