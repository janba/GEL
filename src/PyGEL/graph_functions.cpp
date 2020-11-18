//
//  graph_functions.cpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//


#include <iostream>
#include <string>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/Geometry/graph_util.h>
#include "Graph.h"
#include "Manifold.h"

#include "graph_functions.h"

using namespace std;
using namespace Geometry;
using namespace HMesh;

bool graph_load(Graph_ptr _g_ptr, const char* _file_name) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    const string file_name(_file_name);
    *g_ptr = graph_load(file_name);
    return g_ptr->no_nodes()>0;
}

bool graph_save(Graph_ptr _g_ptr, const char* _file_name) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    const string file_name(_file_name);
    return Geometry::graph_save(file_name, *g_ptr);
}

void graph_to_mesh_cyl(Graph_ptr _g_ptr, Manifold_ptr _m_ptr, float fudge) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    graph_to_mesh_cyl(*g_ptr, *m_ptr, fudge);
}

void graph_smooth(Graph_ptr _g_ptr, const int iter, const float alpha) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    smooth_graph(*g_ptr, iter, alpha);
}

void graph_edge_contract(Graph_ptr _g_ptr, double dist_thresh) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    graph_edge_contract(*g_ptr, dist_thresh);
}

void graph_prune(Graph_ptr _g_ptr) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    prune(*g_ptr);
}

void graph_LS_skeleton(Graph_ptr _g_ptr, Graph_ptr _skel_ptr, bool sampling) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* skel_ptr = reinterpret_cast<AMGraph3D*>(_skel_ptr);
    auto seps = local_separators(*g_ptr, sampling);
    auto [skel, mapping]  = skeleton_from_node_set_vec(*g_ptr, seps);
    *skel_ptr = skel;
}
