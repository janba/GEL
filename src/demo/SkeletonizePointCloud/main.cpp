

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

int main(int argc, char** argv) {
    // Load the points, use a kD tree to connect nearest points.
    // First argument is radius, second is max number of neighbors.
    AMGraph3D g = graph_from_points("./data/armadillo.pts", 7.5, 15);
    
    // Now compute the separators for the skeletonization
    NodeSetVec separators = multiscale_local_separators(g);
    
    // Compute the skeleton from the separators.
    AMGraph3D s = skeleton_from_node_set_vec(g, separators).first;
    
    // And reconstruct a mesh from the skeleton.
    vector<double> radii(s.no_nodes());
    for (auto n: s.node_ids())
        radii[n] = 3;
    Manifold m = graph_to_FEQ(s, radii);
    obj_save("data/armadillo_recon.obj", m);
    return 0;
}
