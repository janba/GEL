//
//  Graph.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/12/2015.
//  Copyright © 2015 J. Andreas Bærentzen. All rights reserved.
//

#include <iostream>
#include <map>
#include <queue>
#include "Graph.h"

namespace Geometry {
    
    using namespace CGLA;
    using namespace std;
    
    void AMGraph3D::merge_nodes(NodeID n0, NodeID n1) {
        CGLA::Vec3d p_new = 0.5*(pos[n0]+pos[n1]);
        
        pos[n0] = Vec3d(CGLA_NAN);
        pos[n1] = p_new;
        
        for(auto n : neighbors(n0)) {
            edge_map[n].erase(n0);
            if(n != n1)
                connect_nodes(n, n1);
        }
        edge_map[n0].clear();
    }

    
    /// Special ID value for invalid node
    const AMGraph3D::NodeID AMGraph::InvalidNodeID = std::numeric_limits<size_t>::max();
    
    /// Special ID value for invalid edge
    const AMGraph3D::EdgeID AMGraph::InvalidEdgeID = std::numeric_limits<size_t>::max();
    
    AMGraph3D clean_graph(const AMGraph3D& g)
    {
        AMGraph3D gn; // new graph
        map<AMGraph::NodeID, AMGraph::NodeID> node_map;
        
        // For all nodes that are not too close to previously visited nodes
        // create a node in the new graph
        for(auto n: g.node_ids())
        {
            if(std::isnan(g.pos[n][0])) {
                node_map[n] = AMGraph::InvalidNodeID;
            } else {
                node_map[n] = gn.add_node(g.pos[n]);
                gn.node_color[node_map[n]] = g.node_color[n];
            }
        }
        
        // For all edges in old graph, create a new edge
        for(auto n: g.node_ids())
            if(node_map[n] != AMGraph::InvalidNodeID)
                for(AMGraph::NodeID& nn: g.neighbors(n)) {
                    AMGraph::EdgeID e = gn.connect_nodes(node_map[n], node_map[nn]);
                    if(gn.valid_edge(e)) {
                        AMGraph::EdgeID e_old = g.find_edge(n, nn);
                        if(g.valid_edge(e_old))
                            gn.edge_color[e] = g.edge_color[e_old];
                        else
                            gn.edge_color[e] = Vec3f(0);
                    }
                }
        
        return gn;
    }
    
    BreadthFirstSearch::BreadthFirstSearch(const AMGraph3D& _g, const Util::AttribVec<AMGraph::NodeID, double>& _dist):
    g_ptr(&_g) {
        pred = Util::AttribVec<AMGraph::NodeID, AMGraph::NodeID>(g_ptr->no_nodes(), AMGraph::InvalidNodeID);
        if(_dist.size() == 0) {
            dist = DistAttribVec(g_ptr->no_nodes(), DBL_MAX);
        }
        else {
            dist = _dist;
            for(auto n: g_ptr->node_ids()) {
                bool is_minimum = true;
                for (auto m: g_ptr->neighbors(n)) {
                    if (dist[m] < dist[n])
                        is_minimum = false;
                }
                if(is_minimum) {
                    pq.push(PrimPQElem(-dist[n], n, AMGraph::InvalidNodeID));
                    front.insert(n);
                }
            }
            
        }
    }
    
    void BreadthFirstSearch::add_init_node(AMGraph::NodeID n, double init_dist) {
        pq.push(PrimPQElem(-init_dist, n, AMGraph::InvalidNodeID));
        dist[n] = init_dist;
        front.insert(n);
    }
    
    bool BreadthFirstSearch::Dijkstra_step() {
        bool did_visit = false;
        while(!pq.empty() && !did_visit) {
            last = pq.top();
            auto n = last.node;
            front.erase(n);
            pq.pop();
            if(last.priority == -dist[n]) {
                visited.insert(n);
                did_visit = true;
                for(auto m: g_ptr->neighbors(n)) {
                    double d = sqrt(g_ptr->sqr_dist(n,m)) - last.priority;
                    if(d < dist[m]) {
                        dist[m] = d;
                        pred[m] = n;
                        pq.push(PrimPQElem(-d, m, n));
                        front.insert(m);
                    }
                }
            }
        }
        return did_visit;
    }
    
    bool BreadthFirstSearch::step() {
        if(!pq.empty()) {
            last = pq.top();
            auto n = last.node;
            front.erase(n);
            pq.pop();
            visited.insert(n);
            for(auto m: g_ptr->neighbors(n))
                if (pred[m]==AMGraph::InvalidNodeID){
                    pred[m] = n;
                    pq.push(PrimPQElem(-dist[m], m, n));
                    front.insert(m);
                }
            return true;
        }
        return false;
    }


    
    AMGraph3D minimum_spanning_tree(const AMGraph3D& g, AMGraph::NodeID root)
    {
        if(root == AMGraph::InvalidNodeID)
            root = 0;
        
        AMGraph3D gn;
        for(auto n: g.node_ids())
            gn.add_node(g.pos[n]);

        BreadthFirstSearch bfs(g);
        bfs.add_init_node(root);
        while(bfs.Dijkstra_step()) {
            auto last = bfs.get_last();
            gn.connect_nodes(last, bfs.get_pred(last));
        }
        
        return gn;
    }
    
 
    std::vector<AMGraph::NodeSet> connected_components(const AMGraph& g, const AMGraph::NodeSet& s) {
        using NodeID = AMGraph::NodeID;
        using NodeSet = AMGraph::NodeSet;
        using NodeQueue = queue<NodeID>;

        NodeSet s_visited;
        vector<AMGraph::NodeSet> component_vec;
        for(auto nf0 : s) {
            if(s_visited.count(nf0)==0)
            {
                s_visited.insert(nf0);
                NodeQueue Q;
                Q.push(nf0);
                NodeSet component;
                while(!Q.empty())
                {
                    NodeID nf = Q.front();
                    Q.pop();
                    component.insert(nf);
                    for(auto nnf: g.neighbors(nf))
                        if(s.count(nnf) >0 && s_visited.count(nnf)==0) {
                            Q.push(nnf);
                            s_visited.insert(nnf);
                        }
                }
                component_vec.push_back(component);
            }
        }
        return component_vec;
    }
    
    double vertex_separator_curvature(const AMGraph3D& g, const AMGraph::NodeSet& separator, const AMGraph::NodeSet& interior) {
        int front_curvature = 0;
        int outside_sum = 0;
        int inside_sum = 0;
//        AMGraph3D::NodeSet inside_set;
//        AMGraph3D::NodeSet outside_set;
        for(auto n: separator) {
            int inside=0;
            int outside=0;
            int in_sep=0;
            for(auto nn: g.neighbors(n)) {
                if(interior.count(nn)) {
                    inside += 1;
//                    inside_set.insert(nn);
                }
                else if(separator.count(nn))
                    in_sep += 1;
                else {
                    outside += 1;
//                    outside_set.insert(nn);
                }
            }
            inside_sum += inside;
            outside_sum += outside;
            int node_curvature = outside-inside;
            front_curvature += sqr(node_curvature);// + sqr(in_sep-2);
        }
//        auto inside_sz = inside_set.size();
//        auto outside_sz = outside_set.size();
//        return (max(inside_sz,outside_sz)/(1e-150 + min(inside_sz,outside_sz)))/separator.size();
        if(inside_sum == 0 || outside_sum == 0)
            return 1e100;
        return static_cast<double>(front_curvature) / separator.size();
    }

    
}
