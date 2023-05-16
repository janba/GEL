//
//  graph_functions.hpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef graph_functions_hpp
#define graph_functions_hpp

#include "IntVector.h"

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

typedef char* Graph_ptr;
typedef char* Manifold_ptr;

#ifdef __cplusplus
extern "C" {
#endif

DLLEXPORT void graph_from_mesh(Manifold_ptr m_ptr, Graph_ptr g_ptr);

DLLEXPORT bool graph_load(Graph_ptr g_ptr, const char* file_name);
DLLEXPORT bool graph_save(Graph_ptr g_ptr, const char* file_name);

DLLEXPORT void graph_to_mesh_cyl(Graph_ptr g_ptr, Manifold_ptr m_ptr, float fudge);

DLLEXPORT void graph_smooth(Graph_ptr g_ptr, const int iter, const float alpha);
DLLEXPORT int graph_edge_contract(Graph_ptr g_ptr, double dist_thresh);
DLLEXPORT void graph_prune(Graph_ptr g_ptr);
DLLEXPORT void graph_saturate(Graph_ptr _g_ptr, int hops, double dist_frac, double rad);

DLLEXPORT void graph_LS_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector_ptr map_ptr, bool sampling=false);

DLLEXPORT void graph_front_skeleton(Graph_ptr g_ptr, Graph_ptr skel_ptr, IntVector_ptr map_ptr, int N_col, double* colors);



#ifdef __cplusplus
}
#endif

#endif /* graph_functions_hpp */
