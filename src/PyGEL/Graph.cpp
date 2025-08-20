//
//  AMGraph3D.cpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//

#include <iostream>
#include <string>
#include <GEL/Geometry/Graph.h>

#include "Graph.h"

using namespace Geometry;
using namespace CGLA;
using namespace std;
using NodeID = Geometry::AMGraph3D::NodeID;

namespace PyGEL {

AMGraph3D* Graph_new(){
    return new AMGraph3D();
}

AMGraph3D* Graph_copy(AMGraph3D* old_self) {
    return new AMGraph3D(*old_self);
}

void Graph_delete(AMGraph3D* self) {
    delete self;
}

std::vector<size_t> Graph_nodes(AMGraph3D* self) {
    std::vector<size_t> nodes;
    nodes.reserve(self->no_nodes());
    for(auto v: self->node_ids())
        nodes.push_back(v);
    return nodes;
}

std::vector<size_t> Graph_neighbors(AMGraph3D* self, size_t n, char mode) {
    std::vector<size_t> nbors;
    const auto& adj_map = self->edges(n);
    nbors.reserve(adj_map.size());
    if (mode=='e')
        for(auto e: adj_map)
            nbors.push_back(size_t(e.second));
    else
        for(auto e: adj_map)
            nbors.push_back(size_t(e.first));
    return nbors;
}

size_t Graph_positions(AMGraph3D* self, double** pos){
    auto N = self->pos.size();
    *pos = reinterpret_cast<double*>(&(self->pos[0]));
    return N;
}

void Graph_clear(AMGraph3D* self){
    self->clear();
}

void Graph_cleanup(AMGraph3D* self){
    self->cleanup();
}

size_t Graph_add_node(AMGraph3D* self, const std::vector<double>& pos){
    return size_t(self->add_node(Vec3d(pos[0],pos[1],pos[2])));
}

void Graph_remove_node(AMGraph3D* self, size_t n){
    self->remove_node(NodeID(n));
}

bool Graph_node_in_use(AMGraph3D* self, size_t n){
    return self->in_use(NodeID(n));
}

size_t Graph_connect_nodes(AMGraph3D* self, size_t n0, size_t n1){
    return size_t(self->connect_nodes(NodeID(n0), NodeID(n1)));
}

void Graph_disconnect_nodes(AMGraph3D* self, size_t n0, size_t n1){
    self->disconnect_nodes(NodeID(n0), NodeID(n1));
}

void Graph_merge_nodes(AMGraph3D* self, size_t n0, size_t n1, bool avg){
    self->merge_nodes(NodeID(n0), NodeID(n1));
}

double Graph_average_edge_length(AMGraph3D* self){
    double r = self->average_edge_length();
    return r;
}

} // namespace PyGEL
