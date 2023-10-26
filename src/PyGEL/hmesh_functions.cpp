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
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include "Graph.h"
#include "Manifold.h"



using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Geometry;

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
                      double error_thresh) {
    quadric_simplify(*(reinterpret_cast<Manifold*>(m_ptr)), keep_fraction, singular_thresh, error_thresh);
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

void volume_preserving_cc_smooth(Manifold_ptr m_ptr, int iter) {
    volume_preserving_cc_smooth(*(reinterpret_cast<Manifold*>(m_ptr)), iter);
}

void regularize_quads(Manifold_ptr m_ptr, float weight, float shrink, int iter) {
    regularize_quads(*(reinterpret_cast<Manifold*>(m_ptr)), weight, shrink, iter);
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

void taubin_smooth(Manifold_ptr m_ptr, int iter) {
    taubin_smooth(*(reinterpret_cast<Manifold*>(m_ptr)), iter);
}

void laplacian_smooth(Manifold_ptr m_ptr, float weight, int iter) {
    laplacian_smooth(*(reinterpret_cast<Manifold*>(m_ptr)), weight, iter);
}

void volumetric_isocontour(Manifold_ptr m_ptr, int x_dim, int y_dim, int z_dim, float* data,
                           double* _pmin, double* _pmax, float tau, bool make_triangles, bool high_is_inside) {
    Vec3i dims(x_dim, y_dim, z_dim);
    const Vec3d pmin = *(reinterpret_cast<Vec3d*>(_pmin));
    const Vec3d pmax = *(reinterpret_cast<Vec3d*>(_pmax));
    XForm xform(pmin, pmax, dims, 0.0);
    RGrid<float> grid(dims);
    memcpy(grid.get(), data, grid.get_size() * sizeof(float));
    volume_polygonize(xform, grid, *(reinterpret_cast<Manifold*>(m_ptr)), tau, make_triangles, high_is_inside);
}


void graph_to_feq(Graph_ptr _g_ptr, Manifold_ptr _m_ptr, double *node_radii, bool symmetrize, bool use_graph_radii) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    vector<double> node_rs;
    const size_t N = g_ptr->no_nodes();
    node_rs.resize(N);

    if (use_graph_radii)
        for(auto n : g_ptr->node_ids()) 
             node_rs[n] = g_ptr->node_color[n][1];
    else 
        for(auto n : g_ptr->node_ids())
            node_rs[n] = node_radii[n];

    *m_ptr = graph_to_FEQ(*g_ptr, node_rs, symmetrize);
}
