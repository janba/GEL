//
//  Graph.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/02/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef Graph_h
#define Graph_h

#include <queue>
#include <map>
#include <set>
#include <unordered_set>
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
    class AMGraph {
    public:
        
        /// ID type for nodes
        using NodeID = size_t;
        
        /// Node Set type
//        using NodeSet = std::unordered_set<NodeID>;
        using NodeSet = std::set<NodeID>;

        /// ID type for edges
        using EdgeID = size_t;
        
        /// The adjacency map class
        using AdjMap = std::map<NodeID, EdgeID>;
        
        /// Special ID value for invalid node
		static const NodeID InvalidNodeID;// = std::numeric_limits<size_t>::max();
        
        /// Special ID value for invalid edge
		static const EdgeID InvalidEdgeID;// = std::numeric_limits<size_t>::max();
        
    protected:
        
        /// The array containing the actual nodes
        std::vector<AdjMap> edge_map;
        
        /// Number of edges.
        size_t no_edges = 0;
        
    public:
        
        /// Clear the graph
        void clear() {
            edge_map.clear();
            no_edges = 0;
        }
        
        /// Return number of nodes
        size_t no_nodes() const {return edge_map.size();}
        
        /// Return whether the node is valid (i.e. in the graph)
        bool valid_node(NodeID n) const { return n < no_nodes();}
        
        /// Return whether an edge is valid (i.e. in the graph)
        bool valid_edge(EdgeID e) const { return e < no_edges;}
        
        /// Returns true if the graph contains no nodes, false otherwise
        bool empty() const {return edge_map.empty();}
        
        /// The range returned can be used in range based for loops over all node ids
        const Util::Range node_ids() const { return Util::Range(0,edge_map.size());}

        /// The range returned can be used in range based for loops over all node ids
        const Util::Range edge_ids() const { return Util::Range(0,no_edges);}

        /// Add a node to the graph
        NodeID add_node() {
            NodeID id = edge_map.size();
            edge_map.push_back(AdjMap());
            return id;
        }
        
        /// Find an edge in the graph given two nodes. Returns InvalidEdgeID if no such edge found.
        EdgeID find_edge(NodeID n0, NodeID n1) const
        {
            if(valid_node(n0) && valid_node(n1)) {
                auto it = edge_map[n0].find(n1);
                if (it != edge_map[n0].end())
                    return it->second;
            }
            return InvalidEdgeID;
        }
        
        /** Connect two nodes. Returns the id of the created edge or InvalidEdgeID if either the nodes
         were not valid or the edge already existed. */
        EdgeID connect_nodes(NodeID n0, NodeID n1)
        {
            if(valid_node(n0) && valid_node(n1) &&
               find_edge(n0, n1) == InvalidEdgeID) {
                size_t id = no_edges++;
                edge_map[n0][n1] = id;
                edge_map[n1][n0] = id;
                return id;
            }
            return InvalidEdgeID;
        }
        
        void disconnect_nodes(NodeID n0, NodeID n1) {
            if(valid_node(n0) && valid_node(n1)) {
                edge_map[n0].erase(n1);
                edge_map[n1].erase(n0);
            }
        }
        
        /// Return the NodeIDs of nodes adjacent to a given node
        std::vector<NodeID> neighbors(NodeID n) const {
            std::vector<NodeID> nbrs(edge_map[n].size());
            int i=0;
            for(auto edge : edge_map[n])
                nbrs[i++] = edge.first;
            return nbrs;
        }
        
        /// Return the edges - map from NodeID to EdgeID of the current node.
        AdjMap edges(NodeID n) const {
            return edge_map[n];
        }
        
        /// Return a vector of shared neighbors.
        std::vector<NodeID> shared_neighbors(AMGraph::NodeID n0, AMGraph::NodeID n1) {
            auto nbrs0 = neighbors(n0);
            auto nbrs1 = neighbors(n1);
            std::vector<AMGraph::NodeID> nbr_isect;
            set_intersection(begin(nbrs0), end(nbrs0), begin(nbrs1), end(nbrs1), back_inserter(nbr_isect));
            return nbr_isect;
        };

        
        /// Return the number of edges incident on a given node.
        size_t valence(NodeID n) const { return edge_map[n].size(); }
        
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
        
        /// Add a node at arbitrary 3D position
        NodeID add_node(const CGLA::Vec3d& p)
        {
            NodeID n = AMGraph::add_node();
            pos[n] = p;
            node_color[n] = CGLA::Vec3f(0);
            return n;
        }
        
        /// Isolate the node and set its position to NAN -- effectively removing it
        void remove_node(NodeID n) {
            // Remove all edges connecting neighboring nodes back to n0
            for(auto& edge: edge_map[n])
                edge_map[edge.first].erase(n);
            // Remove the outgoing aspect of all edges
            edge_map[n].clear();
            pos[n] = CGLA::Vec3d(CGLA::CGLA_NAN);
        }
        
        
        /// Create edge connecting two nodes. Calls the AMGraph::connect_nodes function
        EdgeID connect_nodes(NodeID n0, NodeID n1)
        {
            EdgeID e = AMGraph::connect_nodes(n0,n1);
            if(valid_edge(e))
                edge_color[e] = CGLA::Vec3f(0);
            return e;
        }
     
        /** Merge two nodes, the first is removed and the second inherits all connections,
            the new position becomes the average. */
        void merge_nodes(NodeID n0, NodeID n1, bool avg_pos=true); 
        
        /// Compute sqr distance between two nodes - not necessarily connected.
        double sqr_dist(NodeID n0, NodeID n1) const {
            if(valid_node(n0) && valid_node(n1))
                return CGLA::sqr_length(pos[n0]-pos[n1]);
            else
                return  CGLA::CGLA_NAN;
        }
        
        /// Compute the average edge length
        double average_edge_length() const {
            int i=0;
            double sum_len = 0;
            for(auto n: node_ids())
                for(auto nn: neighbors(n))
                    if(n<nn)
                    {
                        sum_len += sqrt(sqr_dist(n, nn));
                        ++i;
                    }
            return sum_len / i;
        }
        
    };
    
    struct PrimPQElem {
        double priority = 0;
        AMGraph::NodeID node = AMGraph::InvalidNodeID;
        AMGraph::NodeID parent = AMGraph::InvalidNodeID;
        
        PrimPQElem() {}
        PrimPQElem(double _priority, AMGraph::NodeID _node, AMGraph::NodeID _parent):
        priority(_priority), node(_node), parent(_parent) {}
    };
    
    inline bool operator<(const PrimPQElem& p0, const PrimPQElem& p1)
    {
        return p0.priority < p1.priority;
    }
    
    
    class BreadthFirstSearch {
        using DistAttribVec = Util::AttribVec<AMGraph::NodeID, double>;
        const AMGraph3D* g_ptr;
        std::priority_queue<PrimPQElem> pq;
        AMGraph::NodeSet visited, front;
        PrimPQElem last;

        int T;
        
    public:

        DistAttribVec dist;
        Util::AttribVec<AMGraph::NodeID, AMGraph::NodeID> pred;
        Util::AttribVec<AMGraph::NodeID, int> T_in;
        Util::AttribVec<AMGraph::NodeID, int> T_out;

    public:

        BreadthFirstSearch(const AMGraph3D& _g, const DistAttribVec& _dist = DistAttribVec(0));
        
        void add_init_node(AMGraph::NodeID n, double init_dist = 0.0);

        bool Dijkstra_step();
        bool step();
        
//        int get_T_in(AMGraph::NodeID n) const {return T_in[n];}
//        int get_T_out(AMGraph::NodeID n) const {return T_out[n];}
        AMGraph::NodeID get_last() const { return last.node; }
//        AMGraph::NodeID get_pred(AMGraph::NodeID n) const { return pred[n];}
        AMGraph::NodeSet get_front() const { return front; }
        AMGraph::NodeSet get_interior() const { return visited; }
//        double get_dist(AMGraph::NodeID n) const { return dist[n];}
//        const DistAttribVec& get_dist_vec() const { return dist;}
    };
    
    /** Clean up graph, removing unused nodes and edges. */
    AMGraph3D clean_graph(const AMGraph3D& g);
    
    /** Computes the minimum spanning tree of the argument using Prim's algorithm and returns
     the resulting graph. */
    AMGraph3D minimum_spanning_tree(const AMGraph3D&,
                                    AMGraph::NodeID root = AMGraph::InvalidNodeID);

    /** Given a NodeSet s, split s into connected components and return those in a vector */
    std::vector<AMGraph::NodeSet> connected_components(const AMGraph& g, const AMGraph::NodeSet& s);
    
    double vertex_separator_curvature(const AMGraph3D& g, const AMGraph::NodeSet& s,  const Util::AttribVec<AMGraph::NodeID, int>& t_out);
}

#endif /* Graph_h */
