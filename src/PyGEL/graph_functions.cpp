//
//  graph_functions.cpp
//  PyGEL
//
//  Created by Andreas Bærentzen on 16/11/2020.
//  Copyright © 2020 Jakob Andreas Bærentzen. All rights reserved.
//


#include <iostream>
#include <string>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/Geometry/graph_util.h>
#include "Graph.h"
#include "Manifold.h"

#include "graph_functions.h"

using namespace std;
using namespace Geometry;
using namespace HMesh;

namespace PyGEL {

void graph_from_mesh(Manifold_ptr _m_ptr, Graph_ptr _g_ptr) {
    AMGraph3D& g = *reinterpret_cast<AMGraph3D*>(_g_ptr);
    Manifold& m = *reinterpret_cast<Manifold*>(_m_ptr);
    VertexAttributeVector<AMGraph::NodeID> v2n;

    for(auto v : m.vertices())
        v2n[v] = g.add_node(m.pos(v));
    for(auto h: m.halfedges()) {
        Walker w = m.walker(h);
        if(h<w.opp().halfedge())
            g.connect_nodes(v2n[w.opp().vertex()], v2n[w.vertex()]);
    }
}


bool graph_load(Graph_ptr _g_ptr, const std::string& _file_name) {
    AMGraph3D* g_ptr = _g_ptr;
    *g_ptr = Geometry::graph_load(_file_name);
    return g_ptr->no_nodes()>0;
}

bool graph_save(Graph_ptr _g_ptr, const std::string& _file_name) {
    AMGraph3D* g_ptr = _g_ptr;
    return Geometry::graph_save(_file_name, *g_ptr);
}

void graph_to_mesh_cyl(Graph_ptr _g_ptr, Manifold_ptr _m_ptr, float fudge) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    graph_to_mesh_cyl(*g_ptr, *m_ptr, fudge);
}

void graph_to_mesh_iso(Graph_ptr _g_ptr, Manifold_ptr _m_ptr, float fudge, size_t grid_res) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    Manifold* m_ptr = reinterpret_cast<Manifold*>(_m_ptr);
    graph_to_mesh_iso(*g_ptr, *m_ptr, grid_res, fudge, 0.0);
}

void graph_smooth(Graph_ptr _g_ptr, const int iter, const float alpha) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    smooth_graph(*g_ptr, iter, alpha);
}

int graph_edge_contract(Graph_ptr _g_ptr, double dist_thresh) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    return graph_edge_contract(*g_ptr, dist_thresh);
}

void graph_prune(Graph_ptr _g_ptr) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    prune(*g_ptr);
}

void graph_LS_skeleton(Graph_ptr _g_ptr, Graph_ptr _skel_ptr, IntVector_ptr _map_ptr, bool sampling) {
    using IntVector = vector<size_t>;

    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* skel_ptr = reinterpret_cast<AMGraph3D*>(_skel_ptr);
    IntVector* map_ptr = reinterpret_cast<IntVector*>(_map_ptr);
    map_ptr->resize(g_ptr->no_nodes());

    auto seps = local_separators(*g_ptr, sampling);
    auto [skel, mapping]  = skeleton_from_node_set_vec(*g_ptr, seps);
    *skel_ptr = skel;

    for(auto n: g_ptr->node_ids())
        (*map_ptr)[n] = mapping[n];
}

void graph_MSLS_skeleton(Graph_ptr _g_ptr, Graph_ptr _skel_ptr, IntVector_ptr _map_ptr, int grow_thresh) {
    using IntVector = vector<size_t>;

    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* skel_ptr = reinterpret_cast<AMGraph3D*>(_skel_ptr);
    IntVector* map_ptr = reinterpret_cast<IntVector*>(_map_ptr);
    map_ptr->resize(g_ptr->no_nodes());

    auto seps = multiscale_local_separators(*g_ptr, Geometry::SamplingType::Advanced, grow_thresh, 0.1);
    auto [skel, mapping]  = skeleton_from_node_set_vec(*g_ptr, seps);
    *skel_ptr = skel;

    for(auto n: g_ptr->node_ids())
        (*map_ptr)[n] = mapping[n];
}

void graph_saturate(Graph_ptr _g_ptr, int hops, double dist_frac, double rad) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    saturate_graph(*g_ptr, hops, dist_frac, rad);
}

void graph_front_skeleton(Graph_ptr _g_ptr, Graph_ptr _skel_ptr, IntVector_ptr _map_ptr, int N_col, double* colors, int intervals){
    using IntVector = vector<size_t>;
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* skel_ptr = reinterpret_cast<AMGraph3D*>(_skel_ptr);
    IntVector* map_ptr = reinterpret_cast<IntVector*>(_map_ptr);

    const size_t N = g_ptr->no_nodes();
    map_ptr->resize(N);
    
    vector<AttribVecDouble> dvv(N_col);
    for(int i=0;i<N_col; ++i) {
        dvv[i] = AttribVecDouble(g_ptr->no_nodes());
        for(auto n: g_ptr->node_ids())
            dvv[i][n] = colors[N_col*n+i];
    }
    
    auto seps = front_separators(*g_ptr, dvv, intervals);

    auto [skel, mapping]  = skeleton_from_node_set_vec(*g_ptr, seps);
    *skel_ptr = skel;

    for(auto n: g_ptr->node_ids())
        (*map_ptr)[n] = mapping[n];
} 

void graph_combined_skeleton(Graph_ptr _g_ptr, Graph_ptr _skel_ptr, IntVector_ptr _map_ptr, int N_col, double* colors, int intervals){ 
    using IntVector = vector<size_t>;
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* skel_ptr = reinterpret_cast<AMGraph3D*>(_skel_ptr);
    IntVector* map_ptr = reinterpret_cast<IntVector*>(_map_ptr);

    const size_t N = g_ptr->no_nodes();
    map_ptr->resize(N);
    
    vector<AttribVecDouble> dvv(N_col);
    for(int i=0;i<N_col; ++i) {
        dvv[i] = AttribVecDouble(g_ptr->no_nodes());
        for(auto n: g_ptr->node_ids())
            dvv[i][n] = colors[N_col*n+i];
    }
    
    auto seps = combined_separators(*g_ptr, Geometry::SamplingType::Advanced, 256, 0.1, 0, dvv, intervals);

    auto [skel, mapping]  = skeleton_from_node_set_vec(*g_ptr, seps);
    *skel_ptr = skel;

    for(auto n: g_ptr->node_ids())
        (*map_ptr)[n] = mapping[n];
}



void graph_minimum_spanning_tree(Graph_ptr _g_ptr, Graph_ptr _mst_ptr, int _root) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    AMGraph3D* mst_ptr = reinterpret_cast<AMGraph3D*>(_mst_ptr);
    auto root = AMGraph3D::NodeID(_root);
    *mst_ptr = minimum_spanning_tree(*g_ptr, root);
}

void graph_close_chordless_cycles(Graph_ptr _g_ptr, int _root, int hops, double rad) {
    AMGraph3D* g_ptr = reinterpret_cast<AMGraph3D*>(_g_ptr);
    auto root = AMGraph3D::NodeID(_root);
    close_chordless_cycles(*g_ptr, root, hops, rad);
}

} // namespace PyGEL
