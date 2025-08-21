//
//  graph_functions.hpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef graph_functions_hpp
#define graph_functions_hpp

#include <vector>
#include <string>

#include "Graph.h"
#include "Manifold.h"

namespace PyGEL {
    void graph_from_mesh(HMesh::Manifold* m_ptr, Geometry::AMGraph3D* g_ptr);
    bool graph_load(Geometry::AMGraph3D* g_ptr, const std::string& file_name);
    bool graph_save(Geometry::AMGraph3D* g_ptr, const std::string& file_name);

    void graph_to_mesh_cyl(Geometry::AMGraph3D* g_ptr, HMesh::Manifold* m_ptr, float fudge);
    void graph_to_mesh_iso(Geometry::AMGraph3D* g_ptr, HMesh::Manifold* m_ptr, float fudge, size_t grid_res);

    void graph_smooth(Geometry::AMGraph3D* g_ptr, const int iter, const float alpha);
    int graph_edge_contract(Geometry::AMGraph3D* g_ptr, double dist_thresh);
    void graph_prune(Geometry::AMGraph3D* g_ptr);
    void graph_saturate(Geometry::AMGraph3D* g_ptr, int hops, double dist_frac, double rad);

    std::vector<size_t> graph_LS_skeleton(Geometry::AMGraph3D* g_ptr, Geometry::AMGraph3D* skel_ptr, bool sampling=false);
    std::vector<size_t> graph_MSLS_skeleton(Geometry::AMGraph3D* g_ptr, Geometry::AMGraph3D* skel_ptr, int grow_thresh=64);
    std::vector<size_t> graph_front_skeleton(Geometry::AMGraph3D* g_ptr, Geometry::AMGraph3D* skel_ptr, int N_col, const std::vector<double>& colors, int intervals);
    std::vector<size_t> graph_combined_skeleton(Geometry::AMGraph3D* g_ptr, Geometry::AMGraph3D* skel_ptr, int N_col, const std::vector<double>& colors, int intervals);

    void graph_minimum_spanning_tree(Geometry::AMGraph3D* g_ptr, Geometry::AMGraph3D* mst_ptr, int root);
    void graph_close_chordless_cycles(Geometry::AMGraph3D* g_ptr, int root, int hops, double rad);
}

#endif
