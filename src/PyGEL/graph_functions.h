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
#include "IntVector.h"
#include "Graph.h"
#include "Manifold.h"

namespace PyGEL {
    void graph_from_mesh(Manifold_ptr m_ptr, Graph_ptr g_ptr);
    bool graph_load(Graph_ptr g_ptr, const std::string& file_name);
    bool graph_save(Graph_ptr g_ptr, const std::string& file_name);
    
    void graph_to_mesh_cyl(Graph_ptr g_ptr, Manifold_ptr m_ptr, float fudge);
    void graph_to_mesh_iso(Graph_ptr g_ptr, Manifold_ptr m_ptr, float fudge, size_t grid_res);
    
    void graph_smooth(Graph_ptr g_ptr, const int iter, const float alpha);
    int graph_edge_contract(Graph_ptr g_ptr, double dist_thresh);
    void graph_prune(Graph_ptr g_ptr);
    void graph_saturate(Graph_ptr g_ptr, int hops, double dist_frac, double rad);
    
    void graph_LS_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector& map_ptr, bool sampling=false);
    void graph_MSLS_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector& map_ptr, int grow_thresh=64);
    void graph_front_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector& map_ptr, int N_col, const std::vector<double>& colors, int intervals);
    void graph_combined_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector& map_ptr, int N_col, const std::vector<double>& colors, int intervals);
    
    void graph_minimum_spanning_tree(Graph_ptr g_ptr, Graph_ptr mst_ptr, int root);
    void graph_close_chordless_cycles(Graph_ptr g_ptr, int root, int hops, double rad);
}

#endif
