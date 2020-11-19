//
//  graph_abstraction.hpp
//  MeshEditE
//
//  Created by Jakob Andreas Bærentzen on 30/04/2018.
//  Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
//

#ifndef graph_abstraction_hpp
#define graph_abstraction_hpp

#include <vector>
#include <unordered_set>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/Manifold.h>
//#include <GEL/Geometry/bounding_box_tools.h>

namespace Geometry {

    using AttribVecDouble = Util::AttribVec<AMGraph::NodeID, double>;
    using NodeSetUnordered = std::unordered_set<AMGraph::NodeID>;
    using NodeSet = AMGraph::NodeSet;
    using NodeSetVec = std::vector<std::pair<double,NodeSet>>;

    int test_intersection (const AMGraph3D::NodeSet& set1, const AMGraph3D::NodeSet& set2);

    std::vector<NodeSetUnordered> connected_components(const AMGraph& g, const NodeSetUnordered& s);

    AttribVecDouble smooth_dist(const AMGraph3D& g, const AttribVecDouble& _dist, int smooth_iter=0);
    AttribVecDouble projection(const AMGraph3D& g, const CGLA::Vec3d& dir, int smooth_iter=0);
    AttribVecDouble negate_dist(const AMGraph3D& g, const AttribVecDouble& dist_in);
    AttribVecDouble from_extrema(const AMGraph3D& g);

    void saturate_graph(AMGraph3D& g, int hops, double rad);


    void mcm_smooth_graph(AMGraph3D& g, const int iter, const float alpha);
    void smooth_graph(AMGraph3D& g, const int iter, const float alpha);

    int graph_edge_contract(AMGraph3D& g, double dist_thresh);

    void prune(AMGraph3D& g);

    void color_graph_node_sets(AMGraph3D& g, const NodeSetVec& node_set_vec);

    CGLA::Vec3d geometric_median(const std::vector<CGLA::Vec3d>& pts);
    AMGraph3D voxel_graph_from_mesh(HMesh::Manifold& m, int res);

    std::pair<CGLA::Vec3d, double> approximate_bounding_sphere(const AMGraph3D& g, const NodeSetUnordered& s = NodeSetUnordered({}));


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

    std::pair<double,double> graph_H_dist(const AMGraph3D& g0, const AMGraph3D& g1, size_t samples = 10000);
}
#endif /* graph_abstraction_hpp */
