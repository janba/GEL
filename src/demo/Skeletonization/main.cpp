

//
//  main.cpp
//  SkeletonizePointCloud
//
//  This demo program shows how it is possible to create a graph from a collection of points
//  and then skeletonize the graph. Several other types of processing are included for demonstration
//  purposes.
//
//  Created by Andreas Bærentzen on 04/01/2024.
//
#include <random>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/HMesh/skeleton_to_FEQ.h>
#include <GEL/HMesh/HMesh.h>

using namespace Geometry;
using namespace HMesh;
using namespace std;


/** This function filters the radii of the nodes of g using median filtering for the number of iterations given by
 iter. This is sometimes needed if the radii are noisy.*/
void median_filter_node_radii(AMGraph3D& g, int iter, vector<double>& radii) {
    
    for (NodeID node0: g.node_ids()) {
        NodeSet ns = {node0};
        for(int _=0;_<iter;++_) {
            NodeSet ns_new = ns;
            for (NodeID nn: ns) {
                auto nbors = g.neighbors(nn);
                for(NodeID nnn: nbors)
                    ns_new.insert(nnn);
            }
            ns = ns_new;
        }
        vector<double> widths;
        for(NodeID nn: ns) {
            widths.push_back(radii[nn]);
        }
        sort(widths.begin(),widths.end());
        auto N = widths.size();
        if (N%2==1) {
            radii[node0] = widths[N/2];
        }
        else {
            radii[node0] = 0.5*(widths[(N-1)/2]+widths[N/2]);
        }
    }

}


/** This function subdivides the edges of g whose square lengths are greater than thresh. */
void refine_graph_edges(AMGraph3D& g, double thresh) {
    const auto nodes = g.node_ids();
    for(auto n: nodes) {
        for(auto nn: g.neighbors(n)) 
            if (nn<n) {
                if(g.sqr_dist(n, nn) > thresh) {
                    auto p_new = 0.5*(g.pos[n]+g.pos[nn]);
                    auto n_new = g.add_node(p_new);
                    g.disconnect_nodes(n, nn);
                    g.connect_nodes(n, n_new);
                    g.connect_nodes(nn, n_new);
                    
                }
            }
        }
    }


int main(int argc, char** argv) {
    // Load the points, use a kD tree to connect nearest points.
    // First argument is radius, second is max number of neighbors.
    AMGraph3D g = graph_from_points("../../../data/PointClouds/armadillo.pts", 7.5, 15);
    
    // Now compute the separators for the skeletonization
    NodeSetVec separators = multiscale_local_separators(g);
    
    // Compute the skeleton from the separators.
    AMGraph3D s = skeleton_from_node_set_vec(g, separators).first;
    
    // Refine the skeleton by subdividing edges longer than 4.0
    for(int iter=0;iter<3;++iter)
        refine_graph_edges(s, 4.0*4.0);

    // Now, we compute random radii for each node.
    vector<double> radii(s.no_nodes());
    srand(0);
    for (auto n: s.node_ids())
        radii[n] = rand()%7;
    
    // and filter the radii to make them less random.
    median_filter_node_radii(s, 3, radii);
    
    // And reconstruct a mesh from the skeleton.
    Manifold m = graph_to_FEQ(s, radii);
    obj_save("armadillo_recon.obj", m);
    return 0;
}
