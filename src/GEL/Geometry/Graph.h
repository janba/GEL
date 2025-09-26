//
//  Graph.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/02/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef GEOMETRY_GRAPH_H
#define GEOMETRY_GRAPH_H

#include <queue>
#include <vector>
#include <limits>
#include <GEL/CGLA/Vec.h>
#include <GEL/Util/Range.h>
#include <GEL/Util/AttribVec.h>
#include <GEL/Util/AssociativeContainers.h>

namespace Geometry {
    
/** AMGraph means adjacency map graph. It is a simple graph class that is similar to an adjacency list.
 the difference is that the adjacency is given by a map from node id to edge id instead of a list. The
 reason for this is that we can still iterate over all adjacent nodes (the keys) but we can also find
 ids for the concrete edges. This class is not abstract but also not intended for direct usage. There are no
 attributes for nodes or edges and it is not possible to remove nodes or edges from instances of this
 class. Look to the derived AMGraph3D.
 */
    class AMGraph {
    public:
        
        /// ID type for nodes
        using NodeID = size_t;
        
        /// Node Set type
        // FIXME: Ensure we don't rely on ordering and replace this with Util::HashSet
        using NodeSet = Util::BTreeSet<NodeID>;

        /// ID type for edges
        using EdgeID = size_t;
        
        /// The adjacency map class
        using AdjMap = Util::HashMap<NodeID, EdgeID>;

        /// Special ID value for invalid node
		static constexpr NodeID InvalidNodeID = std::numeric_limits<size_t>::max();

        /// Special ID value for invalid edge
		static constexpr EdgeID InvalidEdgeID = std::numeric_limits<size_t>::max();
        
    protected:
        
        /// The array containing the actual nodes
        std::vector<AdjMap> edge_map;
        
        /// Number of edges.
        size_t no_edges_created = 0;
        
    public:
        
        /// Clear the graph
        void clear() {
            edge_map.clear();
            no_edges_created = 0;
        }
        
        /** Return number of nodes. Note that nodes are never removed from edge_map, so the number of nodes is not decremented.
         If you use the derived class AMGraph3D and call cleanup, the number of nodes will  be the actual number. */
        size_t no_nodes() const {return edge_map.size();}
        
        /** Return number of edges created. Note that if nodes were disconnected thereby removing an edge,
         this number is not decremented. If you use the derived class AMGraph3D and call cleanup, the number
         of edges will  be the actual number. */
        size_t no_edges() const {return no_edges_created;}
        
        /// Return whether the node is valid (i.e. in the graph)
        bool valid_node_id(NodeID n) const { return n < no_nodes();}
        
        /// Return whether an edge is valid (i.e. in the graph)
        bool valid_edge_id(EdgeID e) const { return e < no_edges_created;}
        
        /// Returns true if the graph contains no nodes, false otherwise
        bool empty() const {return edge_map.empty();}
        
        /// The range returned can be used in range based for loops over all node ids
        const Util::Range node_ids() const { return Util::Range(0,edge_map.size());}

        /// The range returned can be used in range based for loops over all node ids
        const Util::Range edge_ids() const { return Util::Range(0,no_edges_created);}

        /// Add a node to the graph
        NodeID add_node() {
            const NodeID id = edge_map.size();
            edge_map.emplace_back();
            return id;
        }
        
        /// Find an edge in the graph given two nodes. Returns InvalidEdgeID if no such edge found.
        [[nodiscard]]
        EdgeID find_edge(const NodeID n0, const NodeID n1) const
        {
            if(valid_node_id(n0) && valid_node_id(n1)) {
                if (const auto it = edge_map[n0].find(n1); it != edge_map[n0].end())
                    return it->second;
            }
            return InvalidEdgeID;
        }
        
        /// Connect two nodes. Returns the id of the existing or created edge. Returns InvalidEdgeID if either of
        /// the nodes were not valid.
        EdgeID connect_nodes(const NodeID n0, const NodeID n1)
        {
            if(valid_node_id(n0) && valid_node_id(n1)) {
                size_t id = find_edge(n0, n1);
                if (id == InvalidEdgeID) {
                    id = no_edges_created++;
                    edge_map[n0][n1] = id;
                    edge_map[n1][n0] = id;
                }
                return id;
            }
            return InvalidEdgeID;
        }

        /// Return the NodeIDs of nodes adjacent to a given node lazily
        [[nodiscard]] std::ranges::input_range
        auto neighbors_lazy(const NodeID n) const
        {
            GEL_ASSERT(valid_node_id(n));
            return edge_map[n] | std::views::keys;
        }

        /// Return the NodeIDs of nodes adjacent to a given node
        [[nodiscard]]
        std::vector<NodeID> neighbors(const NodeID n) const {
            auto iter = neighbors_lazy(n);
            return {iter.begin(), iter.end()};
        }
        
        /// Return the edges - map from NodeID to EdgeID of the current node.
        [[nodiscard]]
        std::ranges::input_range auto edges(const NodeID n) const {
            GEL_ASSERT(valid_node_id(n));
            return edge_map[n] | std::views::all;
        }

        [[nodiscard]] std::ranges::input_range
        auto shared_neighbors_lazy(const NodeID n0, const NodeID n1) const
        {
            GEL_ASSERT(valid_node_id(n0) && valid_node_id(n1));
            return edge_map[n0]
                | std::views::keys
                | std::views::filter([this, n1](const auto& e) { return edge_map[n1].contains(e); });
        }
        
        /// Return a vector of shared neighbors.
        [[nodiscard]]
        std::vector<NodeID> shared_neighbors(const NodeID n0, const NodeID n1) const
        {
            auto iter = shared_neighbors_lazy(n0, n1);
            return {iter.begin(), iter.end()};
        }
        
        /// Return the number of edges incident on a given node.
        size_t valence(NodeID n) const { return edge_map[n].size(); }

        void erase_edge(const NodeID n0, const NodeID n1)
        {
            if (valid_node_id(n0) && valid_node_id(n1)) {
                edge_map[n0].erase(n1);
                edge_map[n1].erase(n0);
            }
        }

        void erase_node(const NodeID n0)
        {
            if (valid_node_id(n0)) {
                for (auto neighbor: neighbors_lazy(n0)) {
                    edge_map[neighbor].erase(n0);
                }
                edge_map[n0].clear();
            }
        }
        
    };
    
    
    /** AMGraph3D extends the AMGraph class by providing attributes for position, edge and node
     colors. This is mostly for convenience as these attributes could be stored outside the class,
     but it gives us a concrete 3D graph structure that has many applications. */
    class AMGraph3D final: public AMGraph
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
        
        /// Clean the graph, removing unused nodes and vertices.
        void cleanup();
        
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
        
        /** Returns true if the ID is a valid NodeID and the position is not NaN. These two conditions are true for all nodes
         that are currently in use in the graph. */
        bool in_use(NodeID n) const {
            if (!valid_node_id(n))
                return false;
            if(pos[n] == CGLA::Vec3d(CGLA::CGLA_NAN) && edge_map.empty())
                return false;
            return true;
        }
        
        
        /// Create edge connecting two nodes. Calls the AMGraph::connect_nodes function
        EdgeID connect_nodes(NodeID n0, NodeID n1)
        {
            EdgeID e = AMGraph::connect_nodes(n0,n1);
            if(valid_edge_id(e))
                edge_color[e] = CGLA::Vec3f(0);
            return e;
        }
        
        /** Disconnect nodes. This operation removes the edge from the edge maps of the two formerly connected
         vertices, but the number of edges reported by the super class AMGraph is not decremented, so the edge is only
         invalidated. Call cleanup to finalize removal. */
        void disconnect_nodes(NodeID n0, NodeID n1) {
            if(valid_node_id(n0) && valid_node_id(n1)) {
                edge_map[n0].erase(n1);
                edge_map[n1].erase(n0);
            }
        }
     
        /** Merge two nodes, the first is removed and the second inherits all connections,
            the new position becomes the average. The return value is the node id of the second
         node. */
        NodeID merge_nodes(NodeID n0, NodeID n1, bool avg_pos=true);

        /** Merge all nodes in the vector passed as argument.  The nodes are invalidated, and a
         new node is created at the average position of the nodes given as argument. The NodeID of
         the created node is returned. */
        NodeID merge_nodes(const std::vector<NodeID>& nodes);

        /// Compute sqr distance between two nodes - not necessarily connected.
        double sqr_dist(NodeID n0, NodeID n1) const {
            if(valid_node_id(n0) && valid_node_id(n1))
                return CGLA::sqr_length(pos[n0]-pos[n1]);
            else
                return  CGLA::CGLA_NAN;
        }
        
        /// Compute the average edge length
        double average_edge_length() const {
            unsigned i=0;
            double sum_len = 0;
            for(NodeID n: node_ids())
                for(NodeID nn: neighbors(n))
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
        Util::AttribVec<AMGraph::NodeID, int> mask;
        
    public:

        BreadthFirstSearch(const AMGraph3D& _g, const DistAttribVec& _dist = DistAttribVec(0));
        
        void add_init_node(AMGraph::NodeID n, double init_dist = 0.0);

        bool Dijkstra_step();
        bool step();
        bool Prim_step();

        AMGraph::NodeID get_last() const { return last.node; }
        AMGraph::NodeSet get_front() const { return front; }
        AMGraph::NodeSet get_interior() const { return visited; }
    };
    
    /** Clean up graph, removing unused nodes and edges. */
    AMGraph3D clean_graph(const AMGraph3D& g);
    
    /** Computes the minimum spanning tree of the argument using Prim's algorithm and returns
     the resulting graph. */
    AMGraph3D minimum_spanning_tree(const AMGraph3D&,
                                    AMGraph::NodeID root = AMGraph::InvalidNodeID);

    void close_chordless_cycles(AMGraph3D& g, AMGraph::NodeID root, int hops, double rad);


    /** Given a NodeSet s, split s into connected components and return those in a vector */
    std::vector<AMGraph::NodeSet> connected_components(const AMGraph& g, const AMGraph::NodeSet& s);
    
    double vertex_separator_curvature(const AMGraph3D& g, const AMGraph::NodeSet& s,  const Util::AttribVec<AMGraph::NodeID, int>& t_out);
}

#endif /* Graph_h */
