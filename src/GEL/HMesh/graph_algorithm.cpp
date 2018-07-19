//
//  graph_algorithm.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/05/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#include "graph_algorithm.h"
#include <queue>
#include <utility>

namespace HMesh {

    using namespace std;
    
    DijkstraOutput Dijkstra(const Manifold& m, VertexID v, const VertexSet region)
    {
        DijkstraOutput d_out(m.allocated_vertices());
        
        auto in_region = [&region](VertexID v) {return region.empty() || region.count(v)>0;};

        auto& dist = d_out.dist;
        auto& predecessor = d_out.pred;
        auto& leaves = d_out.leaves;
        
        VertexAttributeVector<int> frozen(m.allocated_vertices(), 0);
        if(in_region(v)) {
            dist[v]=0;
            priority_queue<pair<double,VertexID>> pq;
            pq.push(make_pair(-dist[v], v));
            double max_dist;
            while(!pq.empty())
            {
                VertexID v = pq.top().second;
                max_dist = dist[v];
                pq.pop();
                auto p_v = m.pos(v);
                if(!frozen[v]){
                    frozen[v]=1;
                    
                    bool is_leaf = true;
                    circulate_vertex_ccw(m,v,[&](VertexID vc) {
                        auto p_vc = m.pos(vc);
                        if(in_region(vc) && !frozen[vc]) {
                            double d = dist[v] + length(p_vc - p_v);
                            if(d<dist[vc]) {
                                dist[vc] = d;
                                pq.push(make_pair(-d, vc));
                                predecessor[vc] = v;
                                is_leaf = false;
                            }
                        }
                    });
                    if (is_leaf)
                        leaves.insert(v);
                }
            }
        }
        return d_out;
    }
    
    VertexAttributeVector<int> backpropagate_subtree_sizes(const Manifold& m, const DijkstraOutput& dio) {
        VertexAttributeVector<int> sts(m.allocated_vertices(), 0);
        for (VertexID v0 : dio.leaves) {
            int i=0;
            for (VertexID v = v0; v != InvalidVertexID; v = dio.pred[v],++i)
                sts[v] += i;
        }
        return sts;
        
    }

}
