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
#include <GEL/Geometry/Graph.h>
#include <GEL/Util/AttribVec.h>

namespace Geometry {
    
    using namespace CGLA;
    using namespace std;
    using namespace Util;

    AMGraph3D::NodeID AMGraph3D::merge_nodes(NodeID n0, NodeID n1, bool avg_pos) {
        if (n0 == n1) {
            edge_map[n1].erase(n1);
        }
        else {
            CGLA::Vec3d p_new = pos[n1];
            if(avg_pos) {
                p_new += pos[n0];
                p_new *= 0.5;
            }
            pos[n0] = Vec3d(CGLA_NAN);
            pos[n1] = p_new;
            
            for(auto n:  neighbors(n0)) {
                edge_map[n].erase(n0);
                if(n != n1 && n != n0)
                    connect_nodes(n, n1);
            }
            edge_map[n0].clear();
        }
        return n1;
    }

    AMGraph3D::NodeID AMGraph3D::merge_nodes(const vector<NodeID>& nodes) {
        NodeID n_new = AMGraph::add_node();
        node_color[n_new] = CGLA::Vec3f(0);
        pos[n_new] = Vec3d(0);
        
        NodeSet ns(begin(nodes), end(nodes));
        NodeSet nsn;
        
        for(auto n: ns) {
            pos[n_new] += pos[n];
            for(auto nn : neighbors(n))
                if (ns.count(nn) == 0)
                    nsn.insert(nn);
            edge_map[n].clear();
            pos[n] = Vec3d(CGLA_NAN);
        }
        
        for (auto nn: nsn) {
            connect_nodes(nn, n_new);
            for(auto n: ns)
                edge_map[nn].erase(n);
        }
        
        pos[n_new] /= ns.size();
        return n_new;
    }
    
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
                    if(gn.valid_edge_id(e)) {
                        AMGraph::EdgeID e_old = g.find_edge(n, nn);
                        if(g.valid_edge_id(e_old))
                            gn.edge_color[e] = g.edge_color[e_old];
                        else
                            gn.edge_color[e] = Vec3f(0);
                    }
                }
        
        return gn;
    }

    void AMGraph3D::cleanup() {
        *this = clean_graph(*this);
    }
    
    BreadthFirstSearch::BreadthFirstSearch(const AMGraph3D& _g, const Util::AttribVec<AMGraph::NodeID, double>& _dist):
    g_ptr(&_g) {
        T = 0;
        T_in = Util::AttribVec<AMGraph::NodeID,int>(g_ptr->no_nodes(), INT_MAX);
        T_out = Util::AttribVec<AMGraph::NodeID,int>(g_ptr->no_nodes(), INT_MAX);
        pred = Util::AttribVec<AMGraph::NodeID, AMGraph::NodeID>(g_ptr->no_nodes(), AMGraph::InvalidNodeID);
        mask = Util::AttribVec<AMGraph::NodeID,int>(g_ptr->no_nodes(), 1);
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
//                    cout << dist[n] << " (";
//                    for (auto m: g_ptr->neighbors(n))
//                        cout << dist[m] << " ";
//                    cout << ")" <<endl;
                }
            }
            
        }
    }
    
    void BreadthFirstSearch::add_init_node(AMGraph::NodeID n, double init_dist) {
        pq.push(PrimPQElem(-init_dist, n, AMGraph::InvalidNodeID));
        dist[n] = init_dist;
        front.insert(n);
        T_in[n] = T;
    }
    
    bool BreadthFirstSearch::Dijkstra_step() {
        bool did_visit = false;
        while(!pq.empty() && !did_visit) {
            ++T;
            last = pq.top();
            auto n = last.node;
            front.erase(n);
            T_out[n] = T;
            pq.pop();
            if(last.priority == -dist[n]) {
                visited.insert(n);
                did_visit = true;
                for(auto m: g_ptr->neighbors(n))
                    if(mask[m]) {
                        double d = sqrt(g_ptr->sqr_dist(n,m)) - last.priority;
                        if(d < dist[m]) {
                            dist[m] = d;
                            pred[m] = n;
                            pq.push(PrimPQElem(-d, m, n));
                            front.insert(m);
                            T_in[m] = T;
                    }
                }
            }
        }
        return did_visit;
    }
    
    bool BreadthFirstSearch::step() {
        if(!pq.empty()) {
            auto n = pq.top().node;
            pq.pop();
            if (T_in[n] == INT_MAX) {
                T_in[n] = T;
                front.insert(n);
            }
            T = T + 1;
            T_out[n] = T;
            front.erase(n);
            visited.insert(n);
            for(auto m: g_ptr->neighbors(n))
                if(mask[m])
                    if (T<T_in[m]){
                        pred[m] = n;
                        pq.push(PrimPQElem(-dist[m], m, n));
                        front.insert(m);
                        T_in[m] = T;
                    }
            return true;
        }
        return false;
    }

bool BreadthFirstSearch::Prim_step() {
    bool did_visit = false;
    while(!pq.empty() && !did_visit) {
        ++T;
        last = pq.top();
        auto n = last.node;
        pq.pop();
        if(T < T_out[n]) {
            front.erase(n);
            T_out[n] = T;
            visited.insert(n);
            did_visit = true;
            for(auto m: g_ptr->neighbors(n))
                if(mask[m]) {
                    pred[m] = n;
                    pq.push(PrimPQElem(- g_ptr->sqr_dist(n,m), m, n));
                    front.insert(m);
                    T_in[m] = T;
            }
        }
    }
    return did_visit;
}


    
    AMGraph3D minimum_spanning_tree(const AMGraph3D& g, AMGraph::NodeID root)
    {
        using NodeID = AMGraph3D::NodeID;
        using QElem = tuple<double, NodeID, NodeID>;
        if(root == AMGraph::InvalidNodeID)
            root = 0;
        
        AMGraph3D gn;
        for(auto n: g.node_ids())
            gn.add_node(g.pos[n]);

        AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);

        priority_queue<QElem> Q;
        for (auto n: g.neighbors(root)) {
            auto d = sqr_length(g.pos[n]-g.pos[root]);
            Q.push(make_tuple(-d, root, n));
        }

        while(!Q.empty()) {
            auto [d,n,m] = Q.top();
            Q.pop();
            
            if (!in_tree[m]) {
                in_tree[m] = true;
                gn.connect_nodes(n, m);
                for (auto nn: g.neighbors(m)) {
                    auto d_nn_m = sqr_length(g.pos[nn]-g.pos[m]);
                    Q.push(make_tuple(-d_nn_m, m, nn));
                }
            }
            
        }
            
        return gn;
    }

void close_chordless_cycles(AMGraph3D& g, AMGraph::NodeID root, int hops, double rad)
{
    g.node_color[root] = Vec3f(1,0,0);
    
    double sq_rad = rad*rad;
    
    using NodeID = AMGraph3D::NodeID;
    using QElem = tuple<double, NodeID, NodeID>;
    if(root == AMGraph::InvalidNodeID)
        root = 0;
    
    AttribVec<NodeID, int> cnt(g.no_nodes(), -1);
    AttribVec<NodeID, double> dist(g.no_nodes(), 1e32);
    AttribVec<NodeID, NodeID> parent(g.no_nodes(), AMGraph::InvalidNodeID);
    cnt[root] = 0;
    dist[root] = 0;
    parent[root] = root;
    priority_queue<QElem> Q;
    for (auto n: g.neighbors(root)) {
        auto sq_dist = sqr_length(g.pos[n]-g.pos[root]);
        if (sq_dist < sq_rad) {
            double d = sqrt(sq_dist);
            Q.push(make_tuple(-d, root, n));
            dist[n] = d;
        }
    }
    
    vector<NodeID> bridge_nodes;
    while(!Q.empty()) {
        auto [d,n,m] = Q.top();
        Q.pop();
        
        if (-d == dist[m]) {
            int c = cnt[n]+1;
            parent[m] = n;
            cnt[m] = c;
            bridge_nodes.push_back(m);
            if (c<hops)
                for (auto nn: g.neighbors(m)) {
                    auto sq_dist = sqr_length(g.pos[nn]-g.pos[root]);
                    auto d_nn = length(g.pos[nn]-g.pos[m]) + dist[m];
                    if (sq_dist < sq_rad && d_nn < dist[nn]) {
                        Q.push(make_tuple(-d_nn, m, nn));
                        dist[nn] = d_nn;
                    }
                }
        }
    }
    
    vector<pair<NodeID, NodeID>> loop_pairs;
    for (auto n: bridge_nodes) {
        for (auto m: g.neighbors(n))
            if(parent[m] != AMGraph3D::InvalidNodeID)
                if(n < m && parent[m] != n && parent[n] != m)
                    loop_pairs.push_back(make_pair(n, m));
    }


    auto trace_back = [&](NodeID n) {
        vector<NodeID> path;
        while(n != root) {
            path.push_back(n);
            n = parent[n];
        }
        return path;
    };
//    cout << "Root node: "<< root << endl;
//    cout << "number of loop pairs " << loop_pairs.size() << endl;
    for (auto [n,m]: loop_pairs) {
        auto p0 = trace_back(n);
        auto p1 = trace_back(m);
        if (p0.size()>0 && p1.size()>0 && p0.back()==p1.back())
            continue;
        if (p0.size()+p1.size()+1>3) {
                p0.push_back(root);
                p0.insert(p0.end(), p1.rbegin(), p1.rend());
                bool found_chord = false;
//                cout << "Path len: " << p0.size() << endl;
                for (int i=0; i<p0.size() && !found_chord;++i)
                    for (int j=2; j<p0.size()-1 && !found_chord;++j)
                        {
                            int k=(i+j)%p0.size();
                            if(g.find_edge(p0[i], p0[k]) != AMGraph3D::InvalidEdgeID) {
                                found_chord = true;
//                                cout << "chord " << p0[i] << ", " << p0[k] << " ::: " << i << ", " << j << ", " << k << endl;
//                                for (int ii=0; ii < p0.size() ; ++ii)
//                                    cout << p0[ii] << " ";
//                                cout << endl;
                            }
                        }
//                cout << "Found chord: " << found_chord << endl;
                
                if (!found_chord) {
                    int i = 0;
                    for (int j=2; j<p0.size()-1;++j) {
//                        cout << "Adding chord " << p0.size() << " " << p0[i] << ", " << p0[j] << endl;
                        auto eid = g.connect_nodes(p0[i], p0[j]);
                        g.edge_color[eid] = Vec3f(0,1,0);
                    }
                }
                    
            }
    }
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
                NodeQueue Q;
                Q.push(nf0);
                s_visited.insert(nf0);
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
    
    double vertex_separator_curvature(const AMGraph3D& g, const AMGraph::NodeSet& separator, const Util::AttribVec<AMGraph::NodeID, int>& t_out) {
        int front_curvature = 0;
        int outside_sum = 0;
        int inside_sum = 0;
        AMGraph::NodeSet ns_inside;
        AMGraph::NodeSet ns_outside;

        for(auto n: separator) {
            int inside=0;
            int outside=0;
            int in_sep=0;
            for(auto nn: g.neighbors(n)) {
               if(separator.count(nn))
                    in_sep += 1;
               else if(t_out[nn]< t_out[n]) {
                    inside += 1;
                   ns_inside.insert(nn);
               }
               else {
                    outside += 1;
                   ns_outside.insert(nn);
               }
            }
            inside_sum += inside;
            outside_sum += outside;
            front_curvature += sqr(outside-inside);// + sqr(in_sep-2);
        }
        if(inside_sum == 0 || outside_sum == 0)
            return 1e100;
//        return static_cast<double>(sqr(ns_outside.size()-ns_inside.size())) / separator.size();
        return static_cast<double>(front_curvature) / separator.size();
  }

    
}
