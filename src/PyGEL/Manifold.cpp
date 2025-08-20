//
//  Manifold.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 11/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include <iostream>
#include <string>
#include <GEL/HMesh/HMesh.h>
#include <GEL/HMesh/Delaunay_triangulate.h>

#include "Manifold.h"

using namespace HMesh;
using namespace CGLA;
using namespace std;
namespace py = pybind11;


namespace PyGEL {

Manifold* Manifold_new()
{
    return new Manifold();
}

Manifold* Manifold_from_triangles(const std::vector<double>& vertices, const std::vector<int>& faces) {
    Manifold* m = new Manifold();
    std::vector<int> face_valencies(faces.size() / 3, 3);
    build(*m, vertices.size() / 3, vertices.data(), 
          faces.size() / 3, &face_valencies[0], faces.data());
    return m;
}

Manifold* Manifold_from_points(int N, const std::vector<double>& pts, const Vec& X_axis, const Vec& Y_axis) {
    std::vector<Vec3d> pts3d(N);
    memcpy(pts3d.data(), pts.data(), N * 3 * sizeof(double));
    Manifold* m = new Manifold(Delaunay_triangulate(pts3d, X_axis, Y_axis));
    return m;
}

Manifold* Manifold_copy(Manifold* self)
{
    return new Manifold(*self);
}

void Manifold_merge(Manifold* self, Manifold* other)
{
    self->merge(*other);
}

void Manifold_delete(Manifold* self)
{
    delete self;
}

py::array_t<Scalar> Manifold_positions(Manifold* self) {
    auto N = self->positions_attribute_vector().size();
    py::array_t<Scalar> data(N*3, reinterpret_cast<Scalar*>(self->positions.data()));
    return data;
}

size_t Manifold_no_allocated_vertices(Manifold* self) {
    return self->allocated_vertices();
}

size_t Manifold_no_allocated_faces(Manifold* self) {
    return self->allocated_faces();
}

size_t Manifold_no_allocated_halfedges(Manifold* self) {
    return self->allocated_halfedges();
}

std::vector<size_t> Manifold_vertices(Manifold* self) {
    std::vector<size_t> verts;
    verts.reserve(self->no_vertices());
    for (auto v : self->vertices())
        verts.push_back(v.get_index());
    return verts;
}

std::vector<size_t> Manifold_faces(Manifold* self) {
    std::vector<size_t> faces;
    faces.reserve(self->no_faces());
    for (auto f : self->faces())
        faces.push_back(f.get_index());
    return faces;
}

std::vector<size_t> Manifold_halfedges(Manifold* self) {
    std::vector<size_t> hedges;
    hedges.reserve(self->no_halfedges());
    for (auto h : self->halfedges())
        hedges.push_back(h.get_index());
    return hedges;
}

std::vector<size_t> Manifold_circulate_vertex(Manifold* self, size_t _v, char mode) {
    VertexID v(_v);
    std::vector<size_t> nbrs;
    circulate_vertex_ccw(*self, v, [&](Walker w){
        switch(mode) {
            case 'v': nbrs.push_back(w.vertex().get_index()); break;
            case 'f': nbrs.push_back(w.face().get_index()); break;
            case 'h': nbrs.push_back(w.halfedge().get_index()); break;
        }
    });
    return nbrs;
}

std::vector<size_t> Manifold_circulate_face(Manifold* self, size_t _f, char mode) {
    FaceID f(_f);
    std::vector<size_t> nbrs;
    circulate_face_ccw(*self, f, [&](Walker w){
        switch(mode) {
            case 'v': nbrs.push_back(w.vertex().get_index()); break;
            case 'f': nbrs.push_back(w.opp().face().get_index()); break;
            case 'h': nbrs.push_back(w.halfedge().get_index()); break;
        }
    });
    return nbrs;
}

size_t Manifold_add_face(Manifold* self, const std::vector<double>& pos) {
    auto pts = std::views::iota(0UL, pos.size() / 3) |
        std::views::transform([pos](size_t i) {
            return Vec3d(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
        });
    FaceID f = self->add_face(pts);
    return f.get_index();
}

bool Manifold_remove_face(Manifold* self, size_t fid) {
    return self->remove_face(FaceID(fid));
}

bool Manifold_remove_edge(Manifold* self, size_t hid){
    return self->remove_edge(HalfEdgeID(hid));
}
bool Manifold_remove_vertex(Manifold* self, size_t vid){
    return self->remove_vertex(VertexID(vid));
}

bool Manifold_vertex_in_use(Manifold* self, size_t id){
    return self->in_use(VertexID(id));
}

bool Manifold_face_in_use(Manifold* self, size_t id){
    return self->in_use(FaceID(id));
}

bool Manifold_halfedge_in_use(Manifold* self, size_t id){
    return self->in_use(HalfEdgeID(id));
}

bool Manifold_flip_edge(Manifold* self, size_t _h){
    HalfEdgeID h(_h);
    if(precond_flip_edge(*self, h)) {
        self->flip_edge(h);
        return true;
    }
    return false;
}

bool Manifold_collapse_edge(Manifold* self, size_t _h, bool avg_vertices){
    HalfEdgeID h(_h);
    if(precond_collapse_edge(*self, h)) {
        self->collapse_edge(HalfEdgeID(h),avg_vertices);
        return true;
    }
    return false;
}

size_t Manifold_split_face_by_edge(Manifold* self, size_t f, size_t v0, size_t v1){
    return self->split_face_by_edge(FaceID(f), VertexID(v0), VertexID(v1)).get_index();
}

size_t Manifold_split_face_by_vertex(Manifold* self, size_t f){
    return self->split_face_by_vertex(FaceID(f)).get_index();
}

size_t Manifold_split_edge(Manifold* self, size_t h){
    return self->split_edge(HalfEdgeID(h)).get_index();
}

bool Manifold_stitch_boundary_edges(Manifold* self, size_t h0, size_t h1){
    return self->stitch_boundary_edges(HalfEdgeID(h0), HalfEdgeID(h1));
}
bool Manifold_merge_faces(Manifold* self, size_t f, size_t h){
    return self->merge_faces(FaceID(f),HalfEdgeID(h));
}

size_t Manifold_close_hole(Manifold* self, size_t h){
    return self->close_hole(HalfEdgeID(h)).get_index();
}

void Manifold_cleanup(Manifold* self){
    self->cleanup();
}

// Walker based functions

size_t Walker_next_halfedge(Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.next().halfedge().get_index();
}

size_t Walker_prev_halfedge(Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.prev().halfedge().get_index();
}

size_t Walker_opposite_halfedge(Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.opp().halfedge().get_index();
}
size_t Walker_incident_face(Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.face().get_index();
}
size_t Walker_incident_vertex(Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.vertex().get_index();
}

// Functions we will expose as part of the manifold class interface from Python

bool is_halfedge_at_boundary(const Manifold* m_ptr, size_t _h) {
    return boundary(*m_ptr, HalfEdgeID(_h));
}
bool is_vertex_at_boundary(const Manifold* m_ptr, size_t _v) {
    return boundary(*m_ptr, VertexID(_v));
}

double length(const Manifold* m_ptr, size_t _h) {
    return length(*m_ptr,HalfEdgeID(_h));
}

bool boundary_edge(const Manifold* m_ptr, size_t _v, size_t _h) {
    HalfEdgeID h = boundary_edge(*m_ptr, VertexID(_v));
    _h = h.get_index();
    return h != InvalidHalfEdgeID;
}

size_t valency(const Manifold* m_ptr, size_t _v) {
    return valency(*m_ptr, VertexID(_v));
}

HMesh::Manifold::Vec vertex_normal(const Manifold* m_ptr, size_t _v) {
    return normal(*m_ptr, VertexID(_v));
}

bool connected(const Manifold* m_ptr, size_t _v0, size_t _v1) {
    return connected(*m_ptr, VertexID(_v0), VertexID(_v1));
}

size_t no_edges(const Manifold* m_ptr, size_t _f) {
    return no_edges(*m_ptr,FaceID(_f));
}

HMesh::Manifold::Vec face_normal(const Manifold* m_ptr, size_t _f) {
    return normal(*m_ptr, FaceID(_f));
}

double area(const Manifold* m_ptr, size_t _f) {
    return area(*m_ptr,FaceID(_f));
}

double one_ring_area(const Manifold* m_ptr, size_t _v) {
    return one_ring_area(*m_ptr,VertexID(_v));
}

double mixed_area(const Manifold* m, size_t _v){
    return mixed_area(*m, VertexID(_v));
}

double gaussian_curvature(const Manifold* m, size_t _v){
    return gaussian_curvature(*m, VertexID(_v));
}

double mean_curvature(const Manifold* m, size_t _v){
    return mean_curvature(*m, VertexID(_v));
}

std::vector<double> principal_curvatures(const Manifold* m, size_t _v){
    PrincipalCurvatures pc = principal_curvatures(*m, VertexID(_v));
    std::vector<double> curv_info(8);
    curv_info[0] = pc.min_curvature;
    curv_info[1] = pc.max_curvature;
    curv_info[2] = pc.min_curv_direction[0];
    curv_info[3] = pc.min_curv_direction[1];
    curv_info[4] = pc.min_curv_direction[2];
    curv_info[5] = pc.max_curv_direction[0];
    curv_info[6] = pc.max_curv_direction[1];
    curv_info[7] = pc.max_curv_direction[2];
    return curv_info;
}

double perimeter(const Manifold* m_ptr, size_t _f) {
    return perimeter(*m_ptr,FaceID(_f));
}
HMesh::Manifold::Vec centre(const Manifold* m_ptr, size_t _f) {
    return centre(*m_ptr,FaceID(_f));
}

double total_area(const Manifold* m_ptr) {
    return area(*m_ptr);
}
double volume(const Manifold* m_ptr){
    return volume(*m_ptr);
}

} // namespace PyGEL

namespace PyGEL {
    size_t InvalidIndex = InvalidVertexID.get_index();
}
