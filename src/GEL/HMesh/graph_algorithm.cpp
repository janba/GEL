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
                
                if(!frozen[v]){
                    frozen[v]=1;
                    
                    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw())
                        if(in_region(w.vertex()) && !frozen[w.vertex()])
                        {
                            double d = dist[v] + length(m, w.halfedge());
                            if(d<dist[w.vertex()]) {
                                dist[w.vertex()] = d;
                                pq.push(make_pair(-d, w.vertex()));
                                predecessor[w.vertex()] = v;
                            }
                        }
                }
            }
        }
        return d_out;
    }
    
    

}
