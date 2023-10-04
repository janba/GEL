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
#include "IntVector.h"
#include "Graph.h"

using namespace Geometry;
using namespace CGLA;
using namespace std;
using IntVector = vector<size_t>;
using NodeID = Geometry::AMGraph3D::NodeID;

Graph_ptr Graph_new(){
    return reinterpret_cast<Graph_ptr>(new AMGraph3D());
}

DLLEXPORT Graph_ptr Graph_copy(Graph_ptr _old_self) {
    AMGraph3D* old_self = reinterpret_cast<AMGraph3D*>(_old_self);
    return reinterpret_cast<Graph_ptr>(new AMGraph3D(*old_self));
}


void Graph_delete(Graph_ptr _self) {
    delete reinterpret_cast<AMGraph3D*>(_self);
}

size_t Graph_nodes(Graph_ptr _self, IntVector_ptr _nodes) {
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    IntVector* nodes = reinterpret_cast<IntVector*>(_nodes);
    auto N = self->no_nodes();
    nodes->resize(N);
    size_t i=0;
    for(auto v: self->node_ids())
        (*nodes)[i++] = v;
    return N;
}

size_t Graph_neighbors(Graph_ptr _self, size_t n, IntVector_ptr _nbors, char mode) {
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    IntVector* nbors = reinterpret_cast<IntVector*>(_nbors);
    const auto& adj_map = self->edges(n);
    nbors->resize(adj_map.size());

    int i=0;
    if (mode=='e')
        for(auto e: adj_map)
            (*nbors)[i++] = size_t(e.second);
    else
        for(auto e: adj_map)
            (*nbors)[i++] = size_t(e.first);

    return i;
}



size_t Graph_positions(Graph_ptr _self, double** pos){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    auto N = self->pos.size();
    *pos = reinterpret_cast<double*>(&(self->pos[0]));
    return N;
}

void Graph_clear(Graph_ptr _self){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->clear();
}

void Graph_cleanup(Graph_ptr _self){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->cleanup();
}

size_t Graph_add_node(Graph_ptr _self, double* pos){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    return size_t(self->add_node(Vec3d(pos[0],pos[1],pos[2])));
}

void Graph_remove_node(Graph_ptr _self, size_t n){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->remove_node(NodeID(n));
}

void Graph_node_in_use(Graph_ptr _self, size_t n){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->in_use(NodeID(n));
}

size_t Graph_connect_nodes(Graph_ptr _self, size_t n0, size_t n1){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    return size_t(self->connect_nodes(NodeID(n0), NodeID(n1)));
}

void Graph_disconnect_nodes(Graph_ptr _self, size_t n0, size_t n1){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->disconnect_nodes(NodeID(n0), NodeID(n1));
}

void Graph_merge_nodes(Graph_ptr _self, size_t n0, size_t n1, bool avg){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    self->merge_nodes(NodeID(n0), NodeID(n1));
}

double Graph_average_edge_length(Graph_ptr _self){
    AMGraph3D* self = reinterpret_cast<AMGraph3D*>(_self);
    double r = self->average_edge_length();
    cout << "R (c++) " << r << endl;
    return r;
}
