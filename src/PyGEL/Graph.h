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
    using NodeID = Geometry::AMGraph3D::NodeID;

    Geometry::AMGraph3D* Graph_new();
    Geometry::AMGraph3D* Graph_copy(Geometry::AMGraph3D* self);
    void Graph_delete(Geometry::AMGraph3D* self);
    void Graph_clear(Geometry::AMGraph3D* self);

    std::vector<size_t> Graph_nodes(Geometry::AMGraph3D* self);
    std::vector<size_t> Graph_neighbors(Geometry::AMGraph3D* self, size_t n, char mode='n');

    void Graph_cleanup(Geometry::AMGraph3D* self);
    std::vector<double> Graph_positions(Geometry::AMGraph3D* self);
    double Graph_average_edge_length(Geometry::AMGraph3D* self);

    size_t Graph_add_node(Geometry::AMGraph3D* self, const std::vector<double>& pos);
    void Graph_remove_node(Geometry::AMGraph3D* self, size_t n);
    bool Graph_node_in_use(Geometry::AMGraph3D* self, size_t n);
    size_t Graph_connect_nodes(Geometry::AMGraph3D* self, size_t n0, size_t n1);
    void Graph_disconnect_nodes(Geometry::AMGraph3D* self, size_t n0, size_t n1);
    void Graph_merge_nodes(Geometry::AMGraph3D* self, size_t n0, size_t n1, bool avg);
}


#endif /* Graph_hpp */
