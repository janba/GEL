//
//  graph_io.hpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 25/01/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//

#ifndef graph_io_hpp
#define graph_io_hpp

#include <string>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>

namespace  Geometry {

    /**
     @brief Load a graph from a file.
     @param file_name
     
     Graphs are assumed to be stored in a simple format where:
     - each vertex (aka node) is stored on a single line starting with the character 'n' followed
     by three decimal numbers that together give the 3D spatial coordinates of the node.
     - each edge (aka connection) is stored on a line starting with the character 'c' followed by
     two numbers which are interpreted as the indices of the nodes connected by the edge.
     node numbers are assumed to start from 0.
     */
    AMGraph3D graph_load(const std::string& file_name);

    /**
     @brief Save a graph to a file.
     @param file_name
     
     The graphs are saved to the same format as described above for the graph_load function.
     */
    bool graph_save(const std::string& file_name, const AMGraph3D& g);

    /**
     @brief Load a point set from a file and convert to a graph.
     @param file_name
     @param rad the radius within which we connect two points with an edge
     @param N_closest the maximum number of points which we connect to.

     The points are assumed to be stored in a  very simple text format where each
     point is stored as three numbers on a separate line.*/
    AMGraph3D graph_from_points(const std::string& file_name, double rad, int N_closest);

    /**
     @brief Convert a graph to a mesh using convolution surfaces.
     @param g the input graph
     @param m the output mesh
     @param grid_res the resolution of the voxel grid used in the conversion
     @param fudge is the number added to node size
     @param tau is the threshold for isosurface extraction
     
     This function converts a graph to a skeleton by way of a convolution surface sampled on a voxel grid.
     So where do we get the radius of each node. Well, this is currently stored as the green color in the color associated with each node.
     That is an unfortunate hack, but there you are.
     */
    void graph_to_mesh_iso(const AMGraph3D& g, HMesh::Manifold& m, size_t grid_res, float fudge, float tau);

    /**
    @brief Convert a graph to a mesh using cone stubs.
    @param g the input graph
    @param m the output mesh
    @param fudge is the number added to node size

    This function converts a graph to a skeleton by representing each edge as a cone stub.
    So where do we get the radius of each node. Well, this is currently stored as the green color in the color associated with each node.
    That is an unfortunate hack, but there you are.
    */
    void graph_to_mesh_cyl(const AMGraph3D& g, HMesh::Manifold& m, float fudge);

    /** Save a graph in the ply format. This simplistic ply saver actually does not output the mesh connectivity but just
     the point cloud. The raison d'etre is that some codes for skeletonization of point clouds take ply files as input, so we
     need it for comparison. */
    bool graph_save_ply(const std::string& fn, const AMGraph3D& g);

}
#endif /* graph_io_hpp */
