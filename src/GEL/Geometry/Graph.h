//
//  Graph.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/02/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef Graph_h
#define Graph_h

#include <map>
#include <vector>
#include <limits>
#include "../CGLA/Vec3d.h"
#include "../Util/Range.h"
#include "../Util/AttribVec.h"

namespace Geometry {
    
    /** AMGraph means adjacency map graph. It is a simple graph class that is similar to an adjacency list.
     the difference is that the adjacency is given by a map from node id to edge id instead of a list. The
     reason for this is that we can still iterate over all adjacent nodes (the keys) but we can also find
     ids for the concrete edges. This class can be used but has no properties for the edges or nodes. Using
     an external data structure, such attributes can be added.
     */
    class AMGraph{
    public:
        
        /// ID type for nodes
        using NodeID = size_t;
        
        /// ID type for edges
        using EdgeID = size_t;
        
        /// The adjacency map class
        using AdjMap = std::map<NodeID, EdgeID>;
        
        /// Special ID value for invalid node
        static const NodeID InvalidNodeID = std::numeric_limits<size_t>::max();
        
        /// Special ID value for invalid edge
        static const EdgeID InvalidEdgeID = std::numeric_limits<size_t>::max();
        
    protected:
        
        /// The array containing the actual nodes
        std::vector<AdjMap> nodes;
        
        /// Number of edges.
        size_t no_edges = 0;
        
    public:
        
        /// Clear the graph
        void clear() {
            nodes.clear();
            no_edges = 0;
        }
        
        /// Return number of nodes
        size_t no_nodes() const {return nodes.size();}
        
        /// Return whether the node is valid (i.e. in the graph)
        bool valid_node(NodeID n) const { return n < no_nodes();}
        
        /// Return whether an edge is valid (i.e. in the graph)
        bool valid_edge(EdgeID e) const { return e < no_edges;}
        
        /// Returns true if the graph contains no nodes, false otherwise
        bool empty() const {return nodes.empty();}
        
        /// The range returned can be used in range based for loops over all node ids
        const Util::Range node_ids() const { return Util::Range(0,nodes.size());}
        
        /// Add a node to the graph
        NodeID add_node() {
            NodeID id = nodes.size();
            nodes.push_back(AdjMap());
            return id;
        }
        
        
        /// Find an edge in the graph given two nodes. Returns InvalidEdgeID if no such edge found.
        EdgeID find_edge(NodeID n0, NodeID n1) const
        {
            if(valid_node(n0)) {
                for(auto p: nodes[n0])
                    if(p.first == n1)
                        return p.second;
            }
            return InvalidEdgeID;
        }
        
        /** Connect two nodes. Returns the id of the created edge or InvalidEdgeID if either the nodes
         were not valid or the edge already existed. */
        EdgeID connect_nodes(NodeID n0, NodeID n1)
        {
            if(valid_node(n0) && valid_node(n1) &&
               find_edge(n0, n1) == InvalidNodeID) {
                size_t id = no_edges++;
                nodes[n0][n1] = id;
                nodes[n1][n0] = id;
                return id;
            }
            return InvalidEdgeID;
        }
        
        /// Return the adjacency map for a given node
        AdjMap neighbors(NodeID n) const {return nodes[n];}
        
        /// Return the number of edges incident on a given node.
        size_t valence(NodeID n) const { return nodes[n].size(); }
        
    };
    
    
    /** AMGraph3D extends the AMGraph class by providing attributes for position, edge and node
     colors. This is mostly for convenience as these attributes could be stored outside the class,
     but it gives us a concrete 3D graph structure that has many applications. */
    class AMGraph3D: public AMGraph
    {
    public:
        
        /// position attribute for each node
        Util::AttribVec<AMGraph::NodeID, CGLA::Vec3d> pos;
        
        /// Edge color for each edge
        Util::AttribVec<AMGraph::EdgeID, CGLA::Vec3f> edge_color;
        
        /// Node colors for each node
        Util::AttribVec<AMGraph::NodeID, CGLA::Vec3f> node_color;
        
        /// Clear the graph
        void clear()
        {
            AMGraph::clear();
            pos.clear();
            edge_color.clear();
            node_color.clear();
        }
        
        /// Add a node, initializing position to (0,0,0)
        NodeID add_node()
        {
            return add_node(CGLA::Vec3d(0));
        }
        
        /// Add a node at arbitrary 3D position
        NodeID add_node(const CGLA::Vec3d& p)
        {
            NodeID n = AMGraph::add_node();
            pos[n] = p;
            node_color[n] = CGLA::Vec3f(0);
            return n;
        }
        
        
        /// Create edge connecting two nodes. Calls the AMGraph::connect_nodes function
        EdgeID connect_nodes(NodeID n0, NodeID n1)
        {
            EdgeID e = AMGraph::connect_nodes(n0,n1);
            if(valid_edge(e))
                edge_color[e] = CGLA::Vec3f(0);
            return e;
        }
        
        /// Compute sqr distance between two nodes - not necessarily connected.
        double sqr_dist(NodeID n0, NodeID n1) const {
            if(valid_node(n0) && valid_node(n1))
                return CGLA::sqr_length(pos[n0]-pos[n1]);
            else
                return CGLA::CGLA_NAN;
        }
        
    };
    
    /// Merges all nodes of g within a distance of thresh and returns the resulting graph
    AMGraph3D merge_coincident_nodes(const AMGraph3D& g, double thresh = 1e-12);
    
    /** Computes the minimum spanning tree of the argument using Prim's algorithm and returns
     the resulting graph. */
    AMGraph3D minimum_spanning_tree(const AMGraph3D&);
}
#endif /* Graph_h */
