//
//  Graph.hpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include "IntVector.h"

typedef char* Graph_ptr;

#ifdef __cplusplus
extern "C" {
#endif

DLLEXPORT Graph_ptr Graph_new();
DLLEXPORT Graph_ptr Graph_copy(Graph_ptr self);
DLLEXPORT void Graph_delete(Graph_ptr self);
DLLEXPORT void Graph_clear(Graph_ptr self);

DLLEXPORT size_t Graph_nodes(Graph_ptr self, IntVector_ptr nodes);
DLLEXPORT size_t Graph_neighbors(Graph_ptr _self, size_t n, IntVector_ptr _nbors, char mode='n');


DLLEXPORT void Graph_cleanup(Graph_ptr self);
DLLEXPORT size_t Graph_positions(Graph_ptr self, double** pos);
DLLEXPORT double Graph_average_edge_length(Graph_ptr self);

DLLEXPORT size_t Graph_add_node(Graph_ptr self, double* pos);
DLLEXPORT void Graph_remove_node(Graph_ptr self, size_t n);
DLLEXPORT bool Graph_node_in_use(Graph_ptr self, size_t n);
DLLEXPORT size_t Graph_connect_nodes(Graph_ptr self, size_t n0, size_t n1);
DLLEXPORT void Graph_disconnect_nodes(Graph_ptr self, size_t n0, size_t n1);
DLLEXPORT void Graph_merge_nodes(Graph_ptr self, size_t n0, size_t n1, bool avg);


#ifdef __cplusplus
}
#endif


#endif /* Graph_hpp */
