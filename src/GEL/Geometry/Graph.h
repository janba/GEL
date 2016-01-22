//
//  Graph.hpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/12/2015.
//  Copyright © 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __GEOMETRY_Graph_hpp
#define __GEOMETRY_Graph_hpp

#include <limits>
#include <stdio.h>
#include <vector>

namespace Geometry {
    
    using GraphNodeID = size_t;
    
    GraphNodeID InvalidGraphNodeID = std::numeric_limits<size_t>::max();
    
    using GraphNodeIDVec = std::vector<GraphNodeID>;


    struct GraphNode
    {
        GraphNodeIDVec edges;
        
        GraphNode() {}
        
        bool is_child(GraphNodeID c)
        {
            auto last = std::end(edges);
            if(find(std::begin(edges), last, c) == last)
                return false;
            return true;
        }
        
        void add_edge(GraphNodeID c)
        {
            edges.push_back(c);
        }
    };
    
    template<class T>
    class GraphNodeAttributeVector {
        std::vector<T> items;
        
    public:
         
        T& operator[](GraphNodeID id) {
            if(!(id<items.size()))
                items.resize(id+1);
            return items[id];
        }
        
        const T& operator[](GraphNodeID id) const {
            assert(id<items.size());
            return items[id];
        }
        
    };

    
    class Graph
    {
        bool is_directed;
        std::vector<GraphNode> nodes;
    public:
        size_t no_nodes() const {return nodes.size();}
        GraphNodeID add_node() {
            GraphNode new_node;
            GraphNodeID new_id = nodes.size();
            nodes.push_back(new_node);
            return new_id;
        }
        bool link_directed(GraphNodeID parent_id, GraphNodeID child_id){
            GraphNode& parent = nodes[parent_id];
            GraphNode& child = nodes[child_id];
            if(!parent.is_child(child_id))
               {
                   parent.add_edge(child_id);
                   return true;
               }
            return false;
        }
        
        bool link_undirected(GraphNodeID n0, GraphNodeID n1){
            bool a = link_directed(n0, n1);
            bool b = link_directed(n1, n0);
            return a || b;
        }
    };
    
    
    
    
    
}

#endif /* Graph_hpp */
