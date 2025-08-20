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


namespace PyGEL {
    using namespace Geometry;


    using NodeID = AMGraph3D::NodeID;
    
    AMGraph3D* Graph_new();
    AMGraph3D* Graph_copy(AMGraph3D* self);
    void Graph_delete(AMGraph3D* self);
    void Graph_clear(AMGraph3D* self);

    std::vector<size_t> Graph_nodes(AMGraph3D* self);
    std::vector<size_t> Graph_neighbors(AMGraph3D* self, size_t n, char mode='n');

    void Graph_cleanup(AMGraph3D* self);
    std::vector<double> Graph_positions(AMGraph3D* self);
    double Graph_average_edge_length(AMGraph3D* self);

    size_t Graph_add_node(AMGraph3D* self, const std::vector<double>& pos);
    void Graph_remove_node(AMGraph3D* self, size_t n);
    bool Graph_node_in_use(AMGraph3D* self, size_t n);
    size_t Graph_connect_nodes(AMGraph3D* self, size_t n0, size_t n1);
    void Graph_disconnect_nodes(AMGraph3D* self, size_t n0, size_t n1);
    void Graph_merge_nodes(AMGraph3D* self, size_t n0, size_t n1, bool avg);
}


#endif /* Graph_hpp */
