//
//  Graph.hpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <GEL/Geometry/Graph.h>
#include "IntVector.h"

namespace PyGEL {
    using namespace Geometry;
    using Graph_ptr = AMGraph3D*;
    using IntVector_ptr = IntVector*; // C-style alias
    using NodeID = AMGraph3D::NodeID;
    
    Graph_ptr Graph_new();
    Graph_ptr Graph_copy(Graph_ptr self);
    void Graph_delete(Graph_ptr self);
    void Graph_clear(Graph_ptr self);
    
    size_t Graph_nodes(Graph_ptr self, IntVector& nodes);
    size_t Graph_neighbors(Graph_ptr self, size_t n, IntVector& nbors, char mode='n');
    
    void Graph_cleanup(Graph_ptr self);
    std::vector<double> Graph_positions(Graph_ptr self);
    double Graph_average_edge_length(Graph_ptr self);
    
    size_t Graph_add_node(Graph_ptr self, const std::vector<double>& pos);
    void Graph_remove_node(Graph_ptr self, size_t n);
    bool Graph_node_in_use(Graph_ptr self, size_t n);
    size_t Graph_connect_nodes(Graph_ptr self, size_t n0, size_t n1);
    void Graph_disconnect_nodes(Graph_ptr self, size_t n0, size_t n1);
    void Graph_merge_nodes(Graph_ptr self, size_t n0, size_t n1, bool avg);
}


#endif /* Graph_hpp */
