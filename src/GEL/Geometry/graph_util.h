///
/// @file graph_util.h
/// MeshEditE
///
/// Created by Jakob Andreas Bærentzen on 30/04/2018.
/// Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
///

#ifndef GEOMETRY_GRAPH_UTIL_H
#define GEOMETRY_GRAPH_UTIL_H

#include <vector>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/Manifold.h>
//#include <GEL/Geometry/bounding_box_tools.h>

namespace Geometry {
    using AttribVecDouble = Util::AttribVec<AMGraph::NodeID, double>;
    using NodeSetUnordered = Util::HashSet<AMGraph::NodeID>;
    using NodeSet = AMGraph::NodeSet;
    using NodeSetVec = std::vector<std::pair<double,NodeSet>>;
    using ExpansionMap = std::vector<std::vector<AMGraph::NodeID>>;
    using CapacityVecVec = std::vector<std::vector<size_t>>;
    using NodeQueue = std::queue<AMGraph::NodeID>;

    /// Linear time counting of the number of shared members of set1 and set2.
    int test_intersection (const AMGraph3D::NodeSet& set1, const AMGraph3D::NodeSet& set2);

    /// Returns a vector containing the connected components of set s
    std::vector<NodeSetUnordered> connected_components(const AMGraph& g, const NodeSetUnordered& s);

    /// Returns the connected components of the front of s.
    std::vector<NodeSetUnordered> front_components(const AMGraph3D& g, const NodeSetUnordered & s);

    /// Smooth the attributes in dist associated with graph g smooth_iter times
    AttribVecDouble smooth_dist(const AMGraph3D& g, const AttribVecDouble& dist, int smooth_iter=0);

    /// Project the vertex positions of g onto vector dir and return the attribute vector containing the result smoothed smooth_iter times
    AttribVecDouble projection(const AMGraph3D& g, const CGLA::Vec3d& dir, int smooth_iter=0);

    /// Return the attributes in dist_in but negated.
    AttribVecDouble negate_dist(const AMGraph3D& g, const AttribVecDouble& dist_in);

    /** Add edges to g. For each vertex in g we visit neighbors at a maximum of `hops' graph hops from the original vertex.
     Saturation is here used differently from the conventional graph theoretical usage.. */
    void saturate_graph(AMGraph3D& g, int hops, double dist_frac = 1.0001, double rad = 1e300);

    /// Simple Laplacian graph smoothing. iter specifies number of iterations, and alpha in range [0..1] is the weight.
    void smooth_graph(AMGraph3D& g, const int iter, const float alpha);

    /// Contracts edges in the graph g shorter than dist_thresh. A priority queue is used to contract shorter edges first.
    int graph_edge_contract(AMGraph3D& g, double dist_thresh);

    struct SkeletonPQElem {
        double pri;
        Geometry::AMGraph::NodeID n0, n1;
        SkeletonPQElem(double _pri, Geometry::AMGraph::NodeID _n0, Geometry::AMGraph::NodeID _n1);

        bool operator < (const SkeletonPQElem& pq1) const {
            return pri < pq1.pri;
        }
    };

    /// This function unceremoniously removes leaf vertices if they share an edge with a vertex that has valence greater than two.
    void prune(AMGraph3D& g);

    /// Assign the same (random) color to all nodes in node_set_vec
    void color_graph_node_sets(AMGraph3D& g, const NodeSetVec& node_set_vec);

    /** Computes and returns the point that minimizes the sum of distance to the input points (given in pts).
     A fixed number of iterations is used. */
    CGLA::Vec3d geometric_median(const std::vector<CGLA::Vec3d>& pts);

    /**
     @brief Construct a graph of voxels from input mesh
     @param m is the input mesh
     @param res is the resolution
     @returns a graph whose vertices are interior to m.

        This function computes a distance field from m (at resolution res) and the graph is constructed from the interior voxels.
        each voxel is connected to all 26 neighbors.
     */
     AMGraph3D voxel_graph_from_mesh(HMesh::Manifold& m, int res);

    /** Convert a node set from unordered representation to ordered. This makes it possible to count the number of shared
     nodes in linear time. */
    NodeSet order(const NodeSetUnordered& s);

    /// This function computes the neighbors of s in g. In other words it returns the set of nodes that are connected to s but do not belong to s.
    NodeSetUnordered neighbors(const AMGraph3D& g, const NodeSetUnordered& s);

    /** Compute the approximate bounding sphere for the nodes in graph g optionally restricted to the node set passed as second argument.
        This function returns center and radius of the sphere */
    std::pair<CGLA::Vec3d, double> approximate_bounding_sphere(const AMGraph3D& g, const NodeSetUnordered& s = NodeSetUnordered({}));

    /**
     @brief k means clustering of graph nodes
     @param g is the input graph
     @param N is the desired number of clusters
     @param MAX_ITER is the number of iterations
     @returns a node set vector such that all graph nodes belong to precisely one cluster
     
     This function iteratively clusters vertices of g according to the k-means clustering algorithm.
     each vertex is assigned to the closest cluster (simply using Euclidean distance). To ensure we
     get precisely N clusters, clusters which have a poor boundary to size ratio are ejected when there
     are more than N clusters.
     */
    NodeSetVec k_means_node_clusters(AMGraph3D& g, int N, int MAX_ITER);

    struct LineProj {
        double sqr_dist, t;
    };

    class LineSegment {
        const CGLA::Vec3d p0,p1;
        const CGLA::Vec3d dir;
        const double sqlen;
    public:
        LineSegment(const CGLA::Vec3d& _p0, const CGLA::Vec3d& _p1): p0(_p0), p1(_p1), dir(p1-p0), sqlen(sqr_length(dir)) {}
        
        LineProj sqr_distance(const CGLA::Vec3d& p) const {
            CGLA::Vec3d v =  p - p0;
            double proj = dot(v, dir)/sqlen;
            if(proj < 0)
                return {sqr_length(p-p0),0.0};
            if(proj >1)
                return {sqr_length(p-p1),1.0};
            return { sqr_length(p0 + (p1-p0) * proj - p), proj};
        }
    };

    class GraphDist {
        std::vector<LineSegment> segments;
        double R = 0.0;
        Geometry::KDTree<CGLA::Vec3d, size_t> seg_tree;
        
    public:
        
        GraphDist(const Geometry::AMGraph3D& g);
        double dist(const CGLA::Vec3d& p);
    };

    /** Computes the distance at samples points from graph g0 to g1 and vice versa. H is for Hausdorff. */
    std::pair<double,double> graph_H_dist(const AMGraph3D& g0, const AMGraph3D& g1, size_t samples = 10000);

//    std::pair<int,int> node_symmetry(const AMGraph3D& g,  AMGraph::NodeID n0);
    std::vector<std::pair<int,int>> symmetry_pairs(const AMGraph3D& g, AMGraph::NodeID n, double threshold, bool node_symmetry);
//    void all_symmetry_pairs(AMGraph3D& g, double threshold);
}
#endif /* graph_abstraction_hpp */
