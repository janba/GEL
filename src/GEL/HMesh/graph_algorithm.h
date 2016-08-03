//
//  graph_algorithm.hpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/05/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef HMESH_graph_algorithm_h
#define HMESH_graph_algorithm_h

#include "Manifold.h"

namespace HMesh {
    
    struct DijkstraOutput {
        VertexAttributeVector<double> dist;
        VertexAttributeVector<VertexID> pred;
        
        DijkstraOutput(int n):
        dist(n, DBL_MAX), pred(n, InvalidVertexID) {}
    };
    
    DijkstraOutput Dijkstra(const Manifold& m, VertexID source, VertexSet region = VertexSet());
}


#endif /* graph_algorithm_hpp */
