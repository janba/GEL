#ifndef skeleton_to_FEQ_hpp
#define skeleton_to_FEQ_hpp

#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>

HMesh::Manifold graph_to_FEQ(Geometry::AMGraph3D& g, std::vector<double> node_radii);

#endif
