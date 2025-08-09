//
//  graph_skeletonize.hpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 26/01/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//

#ifndef graph_skeletonize_hpp
#define graph_skeletonize_hpp

#include <GEL/Util/AttribVec.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/CGLA/Vec.h>

namespace Geometry {

    using NodeID = AMGraph::NodeID;
    using NodeSet = AMGraph::NodeSet;
    using NodeSetUnordered = std::unordered_set<NodeID>;
    using NodeSetVec = std::vector<std::pair<double, NodeSet>>;
    using AttribVecDouble = Util::AttribVec<NodeID, double>;
    using ExpansionMap = std::vector<std::vector<AMGraph::NodeID>>;
    using CapacityVecVec = std::vector<std::vector<size_t>>;

    // A set of graphs of different sizes representing the same original graph.
    struct MultiScaleGraph {
        std::vector<AMGraph3D> layers;
        std::vector<ExpansionMap> expansion_map_vec; // expansion_map_vec[layer][nodeID]
        CapacityVecVec capacity_vec_vec; // capacity_vec_vec[layer][nodeID]
    };

    /**
     * @brief Create a multi-scale graph of an input graph.
     * @param g the input graph.
     * @param threshold the size of the smallest layer.
     * @param recursive if the layers are created from the original graph or from the previous layer.
     * @return A multi-scale graph of g.
     *
     * The multi-scale graph is created by simplifying g repeatedly by edge contractions. Each layer is
     * half the size of the previous layer. Using recursive affects how the simplification can be reversed
     * with the expansion map.
     */
    MultiScaleGraph multiscale_graph(const AMGraph3D &g, size_t threshold, bool recursive);

    /**
     @brief Compute separators by marching a front along a scalar field.
     @param g  the graph that we operate on.
     @param dvv  a vector of scalar fields defined on the graph vertices.
     @param intervals an integer which indicates the desired density of separators (lower means higher density)
     @returns This function returns a vector of NodeSets containing a
     number of non-overlapping (local) separators.
     
     Given the input scalar  fields, the function marches along each scalar field and produces
     non-overlapping sets of nodes where each set of nodes represents the front at a
     given point in time.  These are then greedily packed, producing the final
     vector of non-overlapping node sets which is then returned.  This can be seen
     as a mix of Reeb graphs - or as a representation which can be turned into that.
     */
    NodeSetVec front_separators(AMGraph3D &g, const std::vector<AttribVecDouble> &dvv, int intervals);

// SHOULD NOT BE EXPOSED DIRECTLY IN INTERFACE
//    /**
//     For a given graph, g,  and a given node n0, we compute a local separator.
//     The algorithm proceeds in a way similar to Dijkstra, finding a set of nodes separator such that there is anoter set of nodes, front,
//     connected to separator via edges and front consists of two connected components.
//     thick_front indicates whether we want to add a layer of nodes to the front before checking the number of connected  components.
//     persistence is how many iterations the front must have two connected components before we consider the interior
//     a local separator.
//     The growth_threshold sets a maximum size of the separator, with -1 being infinite size.
//     The final node set returned is then thinned to the minimal separator.
//     */
//    Separator
//    local_separator(const AMGraph3D &g, NodeID n0, double quality_noise_level, int optimization_steps,
//                    size_t growth_threshold = -1, const CGLA::Vec3d* static_centre = nullptr);


    enum class SamplingType {
        None, Basic, Advanced
    };

    /**
     @brief Compute a set of local separators from the input graph
     @param g  the graph which is geometrically not modified, but the node colors will
     be changed
     @param quality_noise_level the ratio of smallest to greatest front component size
     needed before we acknowledge that the front has split and a separator been
     found.  The smaller this number, the more spurious branches, but if it is too
     high, we might miss some
     @param optimization_steps indicates the number of times we run a simple
     optimization algorithm based on Dijkstra which aims to make the separator thinner.
     @param sampling chooses the sampling to use. Can either be None, Basic, or Advanced.
     Advanced sampling can sometimes be much faster than Basic but is very dependant on
     picking a good advanced_sampling_threshold.
     @param advanced_sampling_threshold for use together with Advanced sampling.
     Sets the limit on restricted separators when sampling. Higher values might give
     better separators but at the cost of runtime. Optimal value depends on input
     and use case.

     @returns This function returns a vector of NodeSets containing a
     number of non-overlapping (local) separators

     This function finds a number of local separators by using a local front
     propagation method where we iteratively find the point closest to a sphere and
     then expand the sphere to contain the points (when needed).  When the front
     (i.e.  vertices connected to those already found) splits in two, we know that
     at least locally, the found vertices form a separator.  This separator is then
     made minimal by stripping away vertices until no more can be removed without
     rejoining the fronts.  This process is then repeated as many times as required
     by the input parameters, and, finally, the local separators produced are packed
     greedily, and the resulting vector of node sets is returned.
     */
    NodeSetVec local_separators(AMGraph3D &g, SamplingType sampling = SamplingType::None,
                                double quality_noise_level = 0.09,
                                int optimization_steps = 0,
                                size_t advanced_sampling_threshold = 64);

    inline NodeSetVec local_separators(AMGraph3D &g, bool sampling = false,
                                       double quality_noise_level = 0.09,
                                       int optimization_steps = 0) {
        if (sampling) {
            return local_separators(g, SamplingType::Basic, quality_noise_level, optimization_steps);
        } else {
            return local_separators(g, SamplingType::None, quality_noise_level, optimization_steps);
        }
    }

    /**
     @brief Compute a set of local separators from the input graph

     @returns This function returns a vector of NodeSets containing a
        number of non-overlapping (local) separators

     A variations of local_separators the grows restricted separators on a multi-scale graph.
     and transform them into real separators. This method is much faster than local_separators but
     can sometime generate skeletons of slightly lower quality. It also depends on a threshold that
     can be difficult to determine an optimal value for.
     */
    NodeSetVec multiscale_local_separators(AMGraph3D &g, SamplingType sampling = SamplingType::Advanced,
                                size_t grow_threshold = 64,
                                double quality_noise_level = 0.09,
                                int optimization_steps = 0);


    /**
     @brief Convert a vector of (non-overlapping) node sets to a skeleton graph
     @param g the input graph
     @param node_set_vec the set of nodes which we use to produce the skeleton.  Each
     set becomes a single node of the skeleton
     @param merge the resulting skeleton may contain triangles where threee nodes are
     mutually connected forming a clique.  If merge is true, such triangles are
     collapsed into a single vertex joining all the thus connected nodes
     @returns A pair containing the graph which represents the skeleton and a mapping
     from nodes of the original graph to nodes of the skeleton.
     
     Using Dijkstra, we add nodes to each node set in node_set_vec until all nodes are assigned to
     exacly one node set.  The skeleton is then easy to find.  Each  (augmented) node set begets
     a vertex of the skeleton, and whenever two vertices belonging to different node sets are
     connected by an edge, we connect the corresponding skeletal vertices.  Normally, we have
     merge=true and all of the thriangles in the graph are reduced to Steiner-like vertices.  The
     graph structure has a color associated with vertices, and we color the Steiner vertices red.
     The green channel is used to store the estimated radius.  In a better API, we would output these
     things as separate attributes.
     */
    std::pair<AMGraph3D, Util::AttribVec<AMGraph3D::NodeID, AMGraph3D::NodeID>>
    skeleton_from_node_set_vec(AMGraph3D &g,
                               const NodeSetVec &node_set_vec,
                               bool merge = true,
                               int smooth_steps = 0);

    /**
     @brief Convert a set of nodes into a set that partitions the graph
     @param g the input graph
     @param node_set_vec is the vector of input node sets. It is assumed that no node in g
     belongs to more than one of these node sets.
     @returns a node set vector such that all graph nodes belong to precisely one node set
     
     This function simply runs dijkstras algorithm starting from the nodes passed as input.
     Subsequently, each node that did not belong to the input node_set_vec is assigned
     to the node set from which it was reached during the run of Dijkstra
     */
    NodeSetVec maximize_node_set_vec(AMGraph3D &g, const NodeSetVec &node_set_vec);


    /**
     @brief Computes a set of local separators from a graph by combining local separators and front separators
     @returns a node set vector such that all graph nodes belong to precisely one node set
     
     This function runs both the multiscale local separators function and also front separators. Then it combines the
     results.
     */
    NodeSetVec combined_separators(AMGraph3D &g,
                                   SamplingType sampling,
                                   const size_t grow_threshold,
                                   double quality_noise_level,
                                   int optimization_steps,
                                   const std::vector<AttribVecDouble> &dvv,
                                   int intervals);


}
#endif /* graph_skeletonize_hpp */
