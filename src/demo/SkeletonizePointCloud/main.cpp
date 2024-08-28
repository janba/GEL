

//
//  main.cpp
//  SkeletonizePointCloud
//
//  Created by Andreas Bærentzen on 04/01/2024.
//

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/HMesh/skeleton_to_FEQ.h>
#include <GEL/HMesh/HMesh.h>

using namespace Geometry;
using namespace HMesh;
using namespace std;

void median_filter_node_radii(AMGraph3D& g, int iter, vector<double>& radii) {
    
    for (NodeID node0: g.node_ids()) {
        NodeSet ns = {node0};
        for(int _=0;_<iter;++_) {
            for (NodeID nn: ns) {
                auto nbors = g.neighbors(nn);
                if (nbors.size()<3 || nn==node0)    
                    for(NodeID nnn: nbors)
                        if (g.neighbors(nnn).size()<3) // If node is a branch node, we ignore children
                            ns.insert(nnn);
            }
        }
        vector<double> widths;
        for(NodeID nn: ns)
            widths.push_back(radii[nn]);
        sort(widths.begin(),widths.end());
        auto N = widths.size();
        if (N%2==1)
            radii[node0] = widths[N/2];
        else
            radii[node0] = 0.5*(widths[(N-1)/2]+widths[N/2]);
    }
}

int main(int argc, char** argv) {
    // Load the points, use a kD tree to connect nearest points.
    // First argument is radius, second is max number of neighbors.
    AMGraph3D g = graph_from_points("./data/armadillo.pts", 7.5, 15);
    
    // Now compute the separators for the skeletonization
    NodeSetVec separators = multiscale_local_separators(g);
    
    // Compute the skeleton from the separators.
    AMGraph3D s = skeleton_from_node_set_vec(g, separators).first;
    
    // And reconstruct a mesh from the skeleton.
    srand(0);
    vector<double> radii(s.no_nodes());
    for (auto n: s.node_ids())
        radii[n] = s.node_color[n][1];
        // radii[n] *= (double(rand())/RAND_MAX+0.5);
    median_filter_node_radii(s, 2, radii);
    Manifold m = graph_to_FEQ(s, radii);
    obj_save("data/armadillo_recon.obj", m);
    return 0;
}
