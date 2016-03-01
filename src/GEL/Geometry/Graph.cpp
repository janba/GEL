//
//  Graph.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/12/2015.
//  Copyright © 2015 J. Andreas Bærentzen. All rights reserved.
//

#include <map>
#include <queue>
#include "Graph.h"

namespace Geometry {
    
    using namespace CGLA;
    using namespace std;
    
    AMGraph3D merge_coincident_nodes(const AMGraph3D& g, double thresh)
    {
        AMGraph3D gn;
        map<AMGraph::NodeID, AMGraph::NodeID> node_map;
        for(auto n: g.node_ids())
        {
            bool erased = false;
            for(auto m: Util::Range(0,n))
            {
                double d = g.sqr_dist(n,m);
                if(d<thresh) {
                    erased = true;
                    node_map[n] = node_map[m];
                }
            }
            if(!erased)
                node_map[n] = gn.add_node(g.pos[n]);
        }

        for(auto n: g.node_ids())
            for(auto& am: g.neighbors(n))
                gn.connect_nodes(node_map[n], node_map[am.first]);

        return gn;
    }
    
    
    struct PrimPQElem {
        double priority;
        AMGraph::NodeID node;
        AMGraph::NodeID parent;
        
        PrimPQElem(double _priority, AMGraph::NodeID _node, AMGraph::NodeID _parent):
        priority(_priority), node(_node), parent(_parent) {}
    };
    
    bool operator<(const PrimPQElem& p0, const PrimPQElem& p1)
    {
        return p0.priority < p1.priority;
    }

    AMGraph3D minimum_spanning_tree(const AMGraph3D& g)
    {
        
        Vec3d mean(0);
        for(auto n: g.node_ids())
            mean += g.pos[n];
        mean /= g.no_nodes();
        
        double min_d = numeric_limits<double>::max();
        AMGraph::NodeID root = AMGraph::InvalidNodeID;
        for(auto n: g.node_ids())
        {
            double d = sqr_length(g.pos[n] - mean);
            if(d<min_d)
            {
                min_d = d;
                root = n;
            }
        }
        
        AMGraph3D gn;
        
        for(auto n: g.node_ids())
            gn.add_node(g.pos[n]);
        
        priority_queue<PrimPQElem> pq;
        
        for(auto a: g.neighbors(root))
        {
            double d = g.sqr_dist(a.first,root);
            pq.push(PrimPQElem(-d, a.first, root));
        }
        
        vector<bool> visited(g.no_nodes(), false);
        while(!pq.empty())
        {
            auto n = pq.top().node;
            auto p = pq.top().parent;
            pq.pop();
            if(!visited[n]) {
                visited[n] = true;
                gn.connect_nodes(n, p);
                
                for(auto a : g.neighbors(n))
                {
                    AMGraph::NodeID m = a.first;
                    if(!visited[m])
                    {
                        double d = g.sqr_dist(n,m);
                        pq.push(PrimPQElem(-d, m, n));
                    }
                }
            }
        }
        return gn;
    }



}
