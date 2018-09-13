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

Manifold* Manifold_new()
{
    Manifold* m_ptr = new Manifold();
    return m_ptr;
}

Manifold* Manifold_from_triangles(int NV, int NF, double* vertices, int* faces) {
    Manifold* m_ptr = new Manifold();
    vector<int> face_valencies(NF,3);
    m_ptr->build(NV, vertices, NF, &face_valencies[0], faces);
    return m_ptr;
    
}

Manifold* Manifold_from_points(int N, double* pts, double* _X_axis, double* _Y_axis) {
    vector<Vec3d> pts3d(N);
    memcpy(pts3d.data(), pts, N * 3 * sizeof(double));
    Vec3d X_axis(_X_axis[0],_X_axis[1],_X_axis[2]);
    Vec3d Y_axis(_Y_axis[0],_Y_axis[1],_Y_axis[2]);
    Manifold* m_ptr = new Manifold(Delaunay_triangulate(pts3d, X_axis, Y_axis));
    return m_ptr;
}



HMesh::Manifold* Manifold_copy(HMesh::Manifold* self)
{
    Manifold* m2 = new Manifold(*self);
    return m2;
}


void Manifold_delete(Manifold* self)
{
    delete self;
}

size_t Manifold_positions(HMesh::Manifold* self, double** pos) {
    auto N = self->positions_attribute_vector().size();
    *pos = reinterpret_cast<double*>(&(self->positions_attribute_vector().get(VertexID(0))));
    return N;
}

size_t Manifold_no_allocated_vertices(HMesh::Manifold* self) {
    return self->allocated_vertices();
}

size_t Manifold_no_allocated_faces(HMesh::Manifold* self) {
    return self->allocated_faces();

}

size_t Manifold_no_allocated_halfedges(HMesh::Manifold* self) {
    return self->allocated_halfedges();
}


size_t Manifold_vertices(HMesh::Manifold* self, IntVector* verts) {
    auto N = self->no_vertices();
    verts->resize(N);
    size_t i=0;
    for(auto v: self->vertices())
        (*verts)[i++] = v.get_index();
    return N;
}

size_t Manifold_faces(HMesh::Manifold* self, IntVector* faces) {
    auto N = self->no_faces();
    faces->resize(N);
    size_t i=0;
    for(auto f: self->faces())
        (*faces)[i++] = f.get_index();
    return N;
}

size_t Manifold_halfedges(HMesh::Manifold* self, IntVector* hedges) {
    auto N = self->no_halfedges();
    hedges->resize(N);
    size_t i=0;
    for(auto h: self->halfedges())
        (*hedges)[i++] = h.get_index();
    return N;    
}

size_t Manifold_circulate_vertex(HMesh::Manifold* self, size_t _v, char mode, IntVector* nbrs) {
    VertexID v(_v);
    size_t N = circulate_vertex_ccw(*self, v, [&](Walker w){
        switch(mode) {
            case 'v':
                nbrs->push_back(w.vertex().get_index());
                break;
            case 'f':
                nbrs->push_back(w.face().get_index());
                break;
            case 'h':
                nbrs->push_back(w.halfedge().get_index());
                break;
        }
    });
    return N;
}

size_t Manifold_circulate_face(HMesh::Manifold* self, size_t _f, char mode, IntVector* nbrs) {
    FaceID f(_f);
    size_t N = circulate_face_ccw(*self, f, [&](Walker w){
        switch(mode) {
            case 'v':
                nbrs->push_back(w.vertex().get_index());
                break;
            case 'f':
                nbrs->push_back(w.opp().face().get_index());
                break;
            case 'h':
                nbrs->push_back(w.halfedge().get_index());
                break;
        }
    });
    return N;
}

void Manifold_add_face(HMesh::Manifold* self, size_t no_verts, double* pos) {
    vector<Vec3d> pts(no_verts);
    for(size_t i=0;i<no_verts;++i) {
        auto v = Vec3d(pos[3*i],pos[3*i+1],pos[3*i+2]);
        pts[i] = v;
    }
    self->add_face(pts);
}

bool Manifold_remove_face(HMesh::Manifold* self,size_t fid) {
    return self->remove_face(FaceID(fid));
}

bool Manifold_remove_edge(HMesh::Manifold* self,size_t hid){
    return self->remove_edge(HalfEdgeID(hid));
}
bool Manifold_remove_vertex(HMesh::Manifold* self,size_t vid){
    return self->remove_vertex(VertexID(vid));
}

bool Manifold_vertex_in_use(HMesh::Manifold* self,size_t id){
    return self->in_use(VertexID(id));
}

bool Manifold_face_in_use(HMesh::Manifold* self,size_t id){
    return self->in_use(FaceID(id));
}

bool Manifold_halfedge_in_use(HMesh::Manifold* self,size_t id){
    return self->in_use(HalfEdgeID(id));
}

bool Manifold_flip_edge(HMesh::Manifold* self,size_t _h){
    HalfEdgeID h(_h);
    if(precond_flip_edge(*self, h)) {
        self->flip_edge(h);
        return true;
    }
    return false;
}

bool Manifold_collapse_edge(HMesh::Manifold* self,size_t _h, bool avg_vertices){
    HalfEdgeID h(_h);
    if(precond_collapse_edge(*self, h)) {
        self->collapse_edge(HalfEdgeID(h),avg_vertices);
        return true;
    }
    return false;
}

size_t Manifold_split_face_by_edge(HMesh::Manifold* self,size_t f, size_t v0, size_t v1){
    return self->split_face_by_edge(FaceID(f), VertexID(v0), VertexID(v1)).get_index();
}

size_t Manifold_split_face_by_vertex(HMesh::Manifold* self,size_t f){
    return self->split_face_by_vertex(FaceID(f)).get_index();
}

size_t Manifold_split_edge(HMesh::Manifold* self,size_t h){
    return self->split_edge(HalfEdgeID(h)).get_index();
}

bool Manifold_stitch_boundary_edges(HMesh::Manifold* self,size_t h0, size_t h1){
    return self->stitch_boundary_edges(HalfEdgeID(h0), HalfEdgeID(h1));
}
bool Manifold_merge_faces(HMesh::Manifold* self,size_t f, size_t h){
    return self->merge_faces(FaceID(f),HalfEdgeID(h));
}

size_t Manifold_close_hole(HMesh::Manifold* self,size_t h){
    return self->close_hole(HalfEdgeID(h)).get_index();
}

void Manifold_cleanup(HMesh::Manifold* self){
    self->cleanup();
}

// Walker based functions

size_t Walker_next_halfedge(HMesh::Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.next().halfedge().get_index();
}

size_t Walker_prev_halfedge(HMesh::Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.prev().halfedge().get_index();
}


size_t Walker_opposite_halfedge(HMesh::Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.opp().halfedge().get_index();
}
size_t Walker_incident_face(HMesh::Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.face().get_index();
}
size_t Walker_incident_vertex(HMesh::Manifold* m_ptr, size_t _h) {
    Walker w = (*m_ptr).walker(HalfEdgeID(_h));
    return w.vertex().get_index();
}


// Functions we will expose as part of the manifold class size_terface from Python

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
void vertex_normal(const Manifold* m_ptr, size_t _v, CGLA::Vec3d* n) {
    *n = normal(*m_ptr, VertexID(_v));
}
bool connected(const Manifold* m_ptr, size_t _v0, size_t _v1) {
    return connected(*m_ptr, VertexID(_v0), VertexID(_v1));
}

size_t no_edges(const Manifold* m_ptr, size_t _f) {
    return no_edges(*m_ptr,FaceID(_f));
}
void face_normal(const Manifold* m_ptr, size_t _f, CGLA::Vec3d* n) {
    *n = normal(*m_ptr, FaceID(_f));
}
double area(const Manifold* m_ptr, size_t _f) {
    return area(*m_ptr,FaceID(_f));
}
double perimeter(const Manifold* m_ptr, size_t _f) {
    return perimeter(*m_ptr,FaceID(_f));
}
void centre(const Manifold* m_ptr, size_t _f, CGLA::Vec3d* c) {
    *c = centre(*m_ptr,FaceID(_f));
}

