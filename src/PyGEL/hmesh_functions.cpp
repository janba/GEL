//
//  hmesh_functions.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "hmesh_functions.h"
#include <string>

using namespace std;
using namespace HMesh;


bool valid(const Manifold* m_ptr) {
    return valid(*m_ptr);
}
bool closed(const Manifold* m_ptr) {
    return closed(*m_ptr);
}

void bbox(const Manifold* m_ptr, CGLA::Vec3d* pmin, CGLA::Vec3d* pmax) {
    bbox(*m_ptr, *pmin, *pmax);
}
void bsphere(const Manifold* m_ptr, CGLA::Vec3d* c, double* _r) {
    float r;
    bsphere(*m_ptr, *c, r);
    *_r = r;
}

int stitch_mesh(HMesh::Manifold* m_ptr, double rad) {
    return HMesh::stitch_mesh(*m_ptr, rad);
}

bool obj_load(char* fn, HMesh::Manifold* m_ptr) {
    return obj_load(string(fn), *m_ptr, true);
}

bool off_load(char* fn, HMesh::Manifold* m_ptr) {
    return off_load(string(fn), *m_ptr);
}

bool ply_load(char* fn, HMesh::Manifold* m_ptr) {
    return ply_load(string(fn), *m_ptr);
}

bool x3d_load(char* fn, HMesh::Manifold* m_ptr) {
    return x3d_load(string(fn), *m_ptr);
}


bool obj_save(char* fn, HMesh::Manifold* m_ptr) {
    return obj_save(string(fn), *m_ptr);
}

bool off_save(char* fn, HMesh::Manifold* m_ptr) {
    return off_save(string(fn), *m_ptr);

}
bool x3d_save(char* fn, HMesh::Manifold* m_ptr) {
    return x3d_save(string(fn), *m_ptr);
}


void remove_caps(HMesh::Manifold* m_ptr, float thresh) {
    remove_caps(*m_ptr, thresh);
}

void remove_needles(HMesh::Manifold* m_ptr, float thresh, bool averagePositions) {
    remove_needles(*m_ptr, thresh, averagePositions);
}

void close_holes(HMesh::Manifold* m_ptr, int max_size) {
    close_holes(*m_ptr, max_size);
}

void flip_orientation(HMesh::Manifold* m_ptr) {
    flip_orientation(*m_ptr);
}

void merge_coincident_boundary_vertices(HMesh::Manifold* m_ptr, double rad) {
    merge_coincident_boundary_vertices(*m_ptr, rad);
}


void minimize_curvature(HMesh::Manifold* m_ptr, bool anneal) {
    minimize_curvature(*m_ptr, anneal);
}

void minimize_dihedral_angle(HMesh::Manifold* m_ptr, int max_iter, bool anneal, bool alpha, double gamma) {
    minimize_dihedral_angle(*m_ptr, max_iter, anneal, alpha, gamma);
}


void maximize_min_angle(HMesh::Manifold* m_ptr, float thresh, bool anneal) {
    maximize_min_angle(*m_ptr, thresh, anneal);
}

void optimize_valency(HMesh::Manifold* m_ptr, bool anneal) {
    optimize_valency(*m_ptr, anneal);
}

void randomize_mesh(HMesh::Manifold* m_ptr, int max_iter) {
    randomize_mesh(*m_ptr, max_iter);
}

void quadric_simplify(HMesh::Manifold* m_ptr, double keep_fraction,
                      double singular_thresh,
                      bool choose_optimal_positions) {
    quadric_simplify(*m_ptr, keep_fraction, singular_thresh, choose_optimal_positions);
}

float average_edge_length(const HMesh::Manifold* m_ptr) {
    return average_edge_length(*m_ptr);
}

float median_edge_length(const HMesh::Manifold* m_ptr) {
    return median_edge_length(*m_ptr);
}

int refine_edges(HMesh::Manifold* m_ptr, float t) {
    return refine_edges(*m_ptr, t);
}

void cc_split(HMesh::Manifold* m_ptr) {
    cc_split(*m_ptr, *m_ptr);
}

void loop_split(HMesh::Manifold* m_ptr) {
    loop_split(*m_ptr, *m_ptr);
}

void root3_subdivide(HMesh::Manifold* m_ptr) {
    root3_subdivide(*m_ptr, *m_ptr);
}

void rootCC_subdivide(HMesh::Manifold* m_ptr) {
    rootCC_subdivide(*m_ptr, *m_ptr);
}

void butterfly_subdivide(HMesh::Manifold* m_ptr) {
    butterfly_subdivide(*m_ptr, *m_ptr);
}

void cc_smooth(HMesh::Manifold* m_ptr) {
    cc_smooth(*m_ptr);
}

void loop_smooth(HMesh::Manifold* m_ptr) {
    loop_smooth(*m_ptr);
}

void shortest_edge_triangulate(HMesh::Manifold* m_ptr) {
    triangulate(*m_ptr, SHORTEST_EDGE);
}

void ear_clip_triangulate(HMesh::Manifold* m_ptr) {
    triangulate(*m_ptr, CLIP_EAR);
}



