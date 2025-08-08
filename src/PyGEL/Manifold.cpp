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
#include "IntVector.h"
#include "Manifold.h"

using namespace HMesh;
using namespace CGLA;
using namespace std;

namespace PyGEL {

Manifold_ptr Manifold_new()
{
    return new Manifold();
}

Manifold_ptr Manifold_from_triangles(int NV, int NF, double* vertices, int* faces) {
    Manifold* m_ptr = new Manifold();
    vector<int> face_valencies(NF,3);
    build(*m_ptr, NV, vertices, NF, &face_valencies[0], faces);
    return m_ptr;
}

Manifold_ptr Manifold_from_points(int N, double* pts, double* _X_axis, double* _Y_axis) {
    vector<Vec3d> pts3d(N);
    memcpy(pts3d.data(), pts, N * 3 * sizeof(double));
    Vec3d X_axis(_X_axis[0],_X_axis[1],_X_axis[2]);
    Vec3d Y_axis(_Y_axis[0],_Y_axis[1],_Y_axis[2]);
    Manifold* m_ptr = new Manifold(Delaunay_triangulate(pts3d, X_axis, Y_axis));
    return m_ptr;
}

Manifold_ptr Manifold_copy(Manifold_ptr _self)
{
    Manifold* self = _self;
    return new Manifold(*self);
}

void Manifold_merge(Manifold_ptr _self, Manifold_ptr _other)
{
    Manifold* self = _self;
    Manifold* other = reinterpret_cast<Manifold*>(_other);
    self->merge(*other);
}


void Manifold_delete(Manifold_ptr _self)
{
    delete _self;
}

size_t Manifold_positions(Manifold_ptr _self, double** pos) {
    Manifold* self = _self;
    auto N = self->positions_attribute_vector().size();
    *pos = reinterpret_cast<double*>(self->positions.data());
    return N;
}

size_t Manifold_no_allocated_vertices(Manifold_ptr _self) {
    Manifold* self = _self;
    return self->allocated_vertices();
}

size_t Manifold_no_allocated_faces(Manifold_ptr _self) {
    Manifold* self = _self;
    return self->allocated_faces();

}

size_t Manifold_no_allocated_halfedges(Manifold_ptr _self) {
    Manifold* self = _self;
    return self->allocated_halfedges();
}


size_t Manifold_vertices(Manifold_ptr _self, IntVector& verts) {
    Manifold* self = _self;
    auto N = self->no_vertices();
    verts.resize(N);
    size_t i = 0;
    for (auto v : self->vertices())
        verts[i++] = v.get_index();
    return N;
}

size_t Manifold_faces(Manifold_ptr _self, IntVector& faces) {
    Manifold* self = _self;
    auto N = self->no_faces();
    faces.resize(N);
    size_t i = 0;
    for (auto f : self->faces())
        faces[i++] = f.get_index();
    return N;
}

size_t Manifold_halfedges(Manifold_ptr _self, IntVector& hedges) {
    Manifold* self = _self;
    auto N = self->no_halfedges();
    hedges.resize(N);
    size_t i = 0;
    for (auto h : self->halfedges())
        hedges[i++] = h.get_index();
    return N;    
}

size_t Manifold_circulate_vertex(Manifold_ptr _self, size_t _v, char mode, IntVector& nbrs) {
    Manifold* self = _self;
    VertexID v(_v);
    size_t N = circulate_vertex_ccw(*self, v, [&](Walker w){
        switch(mode) {
            case 'v':
                nbrs.push_back(w.vertex().get_index());
                break;
            case 'f':
                nbrs.push_back(w.face().get_index());
                break;
            case 'h':
                nbrs.push_back(w.halfedge().get_index());
                break;
        }
    });
    return N;
}

size_t Manifold_circulate_face(Manifold_ptr _self, size_t _f, char mode, IntVector& nbrs) {
    Manifold* self = _self;
    FaceID f(_f);
    size_t N = circulate_face_ccw(*self, f, [&](Walker w){
        switch(mode) {
            case 'v':
                nbrs.push_back(w.vertex().get_index());
                break;
            case 'f':
                nbrs.push_back(w.opp().face().get_index());
                break;
            case 'h':
                nbrs.push_back(w.halfedge().get_index());
                break;
        }
    });
    return N;
}

size_t Manifold_add_face(Manifold_ptr _self, const std::vector<double>& pos) {
    Manifold* self = _self;

    auto pts = std::views::iota(0UL, pos.size() / 3) |
        std::views::transform([pos](size_t i) {
            return Vec3d(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
        });
    FaceID f = self->add_face(pts);
    return f.get_index();
}

bool Manifold_remove_face(Manifold_ptr _self,size_t fid) {
    Manifold* self = _self;
    return self->remove_face(FaceID(fid));
}

bool Manifold_remove_edge(Manifold_ptr _self,size_t hid){
    Manifold* self = _self;
    return self->remove_edge(HalfEdgeID(hid));
}
bool Manifold_remove_vertex(Manifold_ptr _self,size_t vid){
    Manifold* self = _self;
    return self->remove_vertex(VertexID(vid));
}

bool Manifold_vertex_in_use(Manifold_ptr _self,size_t id){
    Manifold* self = _self;
    return self->in_use(VertexID(id));
}

bool Manifold_face_in_use(Manifold_ptr _self,size_t id){
    Manifold* self = _self;
    return self->in_use(FaceID(id));
}

bool Manifold_halfedge_in_use(Manifold_ptr _self,size_t id){
    Manifold* self = _self;
    return self->in_use(HalfEdgeID(id));
}

bool Manifold_flip_edge(Manifold_ptr _self,size_t _h){
    Manifold* self = _self;
    HalfEdgeID h(_h);
    if(precond_flip_edge(*self, h)) {
        self->flip_edge(h);
        return true;
    }
    return false;
}

bool Manifold_collapse_edge(Manifold_ptr _self,size_t _h, bool avg_vertices){
    Manifold* self = _self;
    HalfEdgeID h(_h);
    if(precond_collapse_edge(*self, h)) {
        self->collapse_edge(HalfEdgeID(h),avg_vertices);
        return true;
    }
    return false;
}

size_t Manifold_split_face_by_edge(Manifold_ptr _self,size_t f, size_t v0, size_t v1){
    Manifold* self = _self;
    return self->split_face_by_edge(FaceID(f), VertexID(v0), VertexID(v1)).get_index();
}

size_t Manifold_split_face_by_vertex(Manifold_ptr _self,size_t f){
    Manifold* self = _self;
    return self->split_face_by_vertex(FaceID(f)).get_index();
}

size_t Manifold_split_edge(Manifold_ptr _self,size_t h){
    Manifold* self = _self;
    return self->split_edge(HalfEdgeID(h)).get_index();
}

bool Manifold_stitch_boundary_edges(Manifold_ptr _self,size_t h0, size_t h1){
    Manifold* self = _self;
    return self->stitch_boundary_edges(HalfEdgeID(h0), HalfEdgeID(h1));
}
bool Manifold_merge_faces(Manifold_ptr _self,size_t f, size_t h){
    Manifold* self = _self;
    return self->merge_faces(FaceID(f),HalfEdgeID(h));
}

size_t Manifold_close_hole(Manifold_ptr _self,size_t h){
    Manifold* self = _self;
    return self->close_hole(HalfEdgeID(h)).get_index();
}

void Manifold_cleanup(Manifold_ptr _self){
    Manifold* self = _self;
    self->cleanup();
}

// Walker based functions

size_t Walker_next_halfedge(Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.next().halfedge().get_index();
}

size_t Walker_prev_halfedge(Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.prev().halfedge().get_index();
}


size_t Walker_opposite_halfedge(Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.opp().halfedge().get_index();
}
size_t Walker_incident_face(Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.face().get_index();
}
size_t Walker_incident_vertex(Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.vertex().get_index();
}


// Functions we will expose as part of the manifold class size_terface from Python

bool is_halfedge_at_boundary(const Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return boundary(*m_ptr, HalfEdgeID(_h));
}
bool is_vertex_at_boundary(const Manifold_ptr _m_ptr, size_t _v) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return boundary(*m_ptr, VertexID(_v));
}

double length(const Manifold_ptr _m_ptr, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return length(*m_ptr,HalfEdgeID(_h));
}

bool boundary_edge(const Manifold_ptr _m_ptr, size_t _v, size_t _h) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    HalfEdgeID h = boundary_edge(*m_ptr, VertexID(_v));
    _h = h.get_index();
    return h != InvalidHalfEdgeID;
}

size_t valency(const Manifold_ptr _m_ptr, size_t _v) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return valency(*m_ptr, VertexID(_v));
}

HMesh::Manifold::Vec vertex_normal(const Manifold_ptr _m_ptr, size_t _v) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return normal(*m_ptr, VertexID(_v));
}

bool connected(const Manifold_ptr _m_ptr, size_t _v0, size_t _v1) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return connected(*m_ptr, VertexID(_v0), VertexID(_v1));
}

size_t no_edges(const Manifold_ptr _m_ptr, size_t _f) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return no_edges(*m_ptr,FaceID(_f));
}

HMesh::Manifold::Vec face_normal(const Manifold_ptr _m_ptr, size_t _f) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return normal(*m_ptr, FaceID(_f));
}

double area(const Manifold_ptr _m_ptr, size_t _f) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return area(*m_ptr,FaceID(_f));
}

double one_ring_area(const Manifold_ptr _m_ptr, size_t _v) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return one_ring_area(*m_ptr,VertexID(_v));
}

double mixed_area(const Manifold_ptr m_ptr, size_t _v){
    Manifold* m = reinterpret_cast<Manifold*>(m_ptr);
    return mixed_area(*m, VertexID(_v));
}

double gaussian_curvature(const Manifold_ptr m_ptr, size_t _v){
    Manifold* m = reinterpret_cast<Manifold*>(m_ptr);
    return gaussian_curvature(*m, VertexID(_v));
}

double mean_curvature(const Manifold_ptr m_ptr, size_t _v){
    Manifold* m = reinterpret_cast<Manifold*>(m_ptr);
    return mean_curvature(*m, VertexID(_v));
}

std::vector<double> principal_curvatures(const Manifold_ptr m_ptr, size_t _v){
    Manifold* m = reinterpret_cast<Manifold*>(m_ptr);
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

double perimeter(const Manifold_ptr _m_ptr, size_t _f) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return perimeter(*m_ptr,FaceID(_f));
}
HMesh::Manifold::Vec centre(const Manifold_ptr _m_ptr, size_t _f) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return centre(*m_ptr,FaceID(_f));
}

double total_area(const Manifold_ptr _m_ptr) {
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return area(*m_ptr);
}
double volume(const Manifold_ptr _m_ptr){
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    return volume(*m_ptr);
}

} // namespace PyGEL

namespace PyGEL {
    size_t InvalidIndex = InvalidVertexID.get_index();
}
