#ifndef skeleton_to_FEQ_hpp
#define skeleton_to_FEQ_hpp

#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>

HMesh::Manifold graph_to_FEQ(const Geometry::AMGraph3D& g, const std::vector<double>& node_radii, bool use_symmetry=true);

void non_rigid_registration(HMesh::Manifold& m, const HMesh::Manifold& m_ref);

#endif
