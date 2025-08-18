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

Graph_ptr Graph_new(){
    return new AMGraph3D();
}

Graph_ptr Graph_copy(Graph_ptr _old_self) {
    AMGraph3D* old_self = _old_self;
    return new AMGraph3D(*old_self);
}

void Graph_delete(Graph_ptr _self) {
    delete _self;
}

std::vector<size_t> Graph_nodes(Graph_ptr _self) {
    AMGraph3D* self = _self;
    std::vector<size_t> nodes;
    nodes.reserve(self->no_nodes());
    for(auto v: self->node_ids())
        nodes.push_back(v);
    return nodes;
}

std::vector<size_t> Graph_neighbors(Graph_ptr _self, size_t n, char mode) {
    AMGraph3D* self = _self;
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

size_t Graph_positions(Graph_ptr _self, double** pos){
    AMGraph3D* self = _self;
    auto N = self->pos.size();
    *pos = reinterpret_cast<double*>(&(self->pos[0]));
    return N;
}

void Graph_clear(Graph_ptr _self){
    AMGraph3D* self = _self;
    self->clear();
}

void Graph_cleanup(Graph_ptr _self){
    AMGraph3D* self = _self;
    self->cleanup();
}

size_t Graph_add_node(Graph_ptr _self, const std::vector<double>& pos){
    AMGraph3D* self = _self;
    return size_t(self->add_node(Vec3d(pos[0],pos[1],pos[2])));
}

void Graph_remove_node(Graph_ptr _self, size_t n){
    AMGraph3D* self = _self;
    self->remove_node(NodeID(n));
}

bool Graph_node_in_use(Graph_ptr _self, size_t n){
    AMGraph3D* self = _self;
    return self->in_use(NodeID(n));
}

size_t Graph_connect_nodes(Graph_ptr _self, size_t n0, size_t n1){
    AMGraph3D* self = _self;
    return size_t(self->connect_nodes(NodeID(n0), NodeID(n1)));
}

void Graph_disconnect_nodes(Graph_ptr _self, size_t n0, size_t n1){
    AMGraph3D* self = _self;
    self->disconnect_nodes(NodeID(n0), NodeID(n1));
}

void Graph_merge_nodes(Graph_ptr _self, size_t n0, size_t n1, bool avg){
    AMGraph3D* self = _self;
    self->merge_nodes(NodeID(n0), NodeID(n1));
}

double Graph_average_edge_length(Graph_ptr _self){
    AMGraph3D* self = _self;
    double r = self->average_edge_length();
    return r;
}

} // namespace PyGEL
