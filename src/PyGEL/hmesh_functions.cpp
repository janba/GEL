//
//  hmesh_functions.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "hmesh_functions.h"
#include <string>
#include <GEL/HMesh/HMesh.h>

using namespace std;
using namespace HMesh;
using namespace CGLA;

bool valid(const Manifold_ptr m_ptr) {
    return valid(*(reinterpret_cast<Manifold*>(m_ptr)));
}
bool closed(const Manifold_ptr m_ptr) {
    return closed(*(reinterpret_cast<Manifold*>(m_ptr)));
}

void bbox(const Manifold_ptr m_ptr, double* pmin, double* pmax) {
    bbox(*(reinterpret_cast<Manifold*>(m_ptr)),
         *(reinterpret_cast<Vec3d*>(pmin)),
         *(reinterpret_cast<Vec3d*>(pmax)));
}
void bsphere(const Manifold_ptr m_ptr, double* c, double* _r) {
    float r;
    bsphere(*(reinterpret_cast<Manifold*>(m_ptr)), *reinterpret_cast<Vec3d*>(c), r);
    *_r = r;
}

int stitch_mesh(Manifold_ptr m_ptr, double rad) {
    return stitch_mesh(*(reinterpret_cast<Manifold*>(m_ptr)), rad);
}

bool load(const char* fn, Manifold_ptr m_ptr) {
    return load(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}


bool obj_load(const char* fn, Manifold_ptr m_ptr) {
    return obj_load(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}

bool off_load(const char* fn, Manifold_ptr m_ptr) {
    return off_load(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}

bool ply_load(const char* fn, Manifold_ptr m_ptr) {
    return ply_load(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}

bool x3d_load(const char* fn, Manifold_ptr m_ptr) {
    return x3d_load(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}


bool obj_save(const char* fn, Manifold_ptr m_ptr) {
    return obj_save(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}

bool off_save(const char* fn, Manifold_ptr m_ptr) {
    return off_save(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));

}
bool x3d_save(const char* fn, Manifold_ptr m_ptr) {
    return x3d_save(string(fn), *(reinterpret_cast<Manifold*>(m_ptr)));
}


void remove_caps(Manifold_ptr m_ptr, float thresh) {
    remove_caps(*(reinterpret_cast<Manifold*>(m_ptr)), thresh);
}

void remove_needles(Manifold_ptr m_ptr, float thresh, bool averagePositions) {
    remove_needles(*(reinterpret_cast<Manifold*>(m_ptr)), thresh, averagePositions);
}

void close_holes(Manifold_ptr m_ptr, int max_size) {
    close_holes(*(reinterpret_cast<Manifold*>(m_ptr)), max_size);
}

void flip_orientation(Manifold_ptr m_ptr) {
    flip_orientation(*(reinterpret_cast<Manifold*>(m_ptr)));
}

void merge_coincident_boundary_vertices(Manifold_ptr m_ptr, double rad) {
    merge_coincident_boundary_vertices(*(reinterpret_cast<Manifold*>(m_ptr)), rad);
}


void minimize_curvature(Manifold_ptr m_ptr, bool anneal) {
    minimize_curvature(*(reinterpret_cast<Manifold*>(m_ptr)), anneal);
}

void minimize_dihedral_angle(Manifold_ptr m_ptr, int max_iter, bool anneal, bool alpha, double gamma) {
    minimize_dihedral_angle(*(reinterpret_cast<Manifold*>(m_ptr)), max_iter, anneal, alpha, gamma);
}


void maximize_min_angle(Manifold_ptr m_ptr, float thresh, bool anneal) {
    maximize_min_angle(*(reinterpret_cast<Manifold*>(m_ptr)), thresh, anneal);
}

void optimize_valency(Manifold_ptr m_ptr, bool anneal) {
    optimize_valency(*(reinterpret_cast<Manifold*>(m_ptr)), anneal);
}

void randomize_mesh(Manifold_ptr m_ptr, int max_iter) {
    randomize_mesh(*(reinterpret_cast<Manifold*>(m_ptr)), max_iter);
}

void quadric_simplify(Manifold_ptr m_ptr, double keep_fraction,
                      double singular_thresh,
                      bool choose_optimal_positions) {
    quadric_simplify(*(reinterpret_cast<Manifold*>(m_ptr)), keep_fraction, singular_thresh, choose_optimal_positions);
}

float average_edge_length(const Manifold_ptr m_ptr) {
    return average_edge_length(*(reinterpret_cast<Manifold*>(m_ptr)));
}

float median_edge_length(const Manifold_ptr m_ptr) {
    return median_edge_length(*(reinterpret_cast<Manifold*>(m_ptr)));
}

int refine_edges(Manifold_ptr m_ptr, float t) {
    return refine_edges(*(reinterpret_cast<Manifold*>(m_ptr)), t);
}

void cc_split(Manifold_ptr m_ptr) {
    cc_split(*(reinterpret_cast<Manifold*>(m_ptr)), *(reinterpret_cast<Manifold*>(m_ptr)));
}

void loop_split(Manifold_ptr m_ptr) {
    loop_split(*(reinterpret_cast<Manifold*>(m_ptr)), *(reinterpret_cast<Manifold*>(m_ptr)));
}

void root3_subdivide(Manifold_ptr m_ptr) {
    root3_subdivide(*(reinterpret_cast<Manifold*>(m_ptr)), *(reinterpret_cast<Manifold*>(m_ptr)));
}

void rootCC_subdivide(Manifold_ptr m_ptr) {
    rootCC_subdivide(*(reinterpret_cast<Manifold*>(m_ptr)), *(reinterpret_cast<Manifold*>(m_ptr)));
}

void butterfly_subdivide(Manifold_ptr m_ptr) {
    butterfly_subdivide(*(reinterpret_cast<Manifold*>(m_ptr)), *(reinterpret_cast<Manifold*>(m_ptr)));
}

void cc_smooth(Manifold_ptr m_ptr) {
    cc_smooth(*(reinterpret_cast<Manifold*>(m_ptr)));
}

void loop_smooth(Manifold_ptr m_ptr) {
    loop_smooth(*(reinterpret_cast<Manifold*>(m_ptr)));
}

void shortest_edge_triangulate(Manifold_ptr m_ptr) {
    triangulate(*(reinterpret_cast<Manifold*>(m_ptr)), SHORTEST_EDGE);
}

void ear_clip_triangulate(Manifold_ptr m_ptr) {
    triangulate(*(reinterpret_cast<Manifold*>(m_ptr)), CLIP_EAR);
}



