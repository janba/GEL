//
//  gem.hpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 24/11/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef face_loop_hmesh
#define face_loop_hmesh

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/Graph.h>

/** FaceLoop data structure. The face loop is represented in terms of halfedges which separate the faces of the loop.
 The data structure contains a lot of descriptive members which are needed for some of the functions that use it. */
struct FaceLoop {
    std::vector<HMesh::HalfEdgeID> hvec;
    HMesh::FaceSet double_cross_faces;
    HMesh::FaceSet interior;
    int interior_faces = 0;
    double area = 0;
    double avg_len = 0;
    double interior_area = 0;
    double cylindricity = 0;
    double valency_imbalance = 0;
    double integral_geodesic_curvature = 0;
    CGLA::Vec3d avg_edge = CGLA::Vec3d(0,0,0);
    CGLA::Vec3d center = CGLA::Vec3d(0,0,0);
    int id=-1;
};

/** Create a hexahedron by adding six faces to the Manifold m passed as input. The other arguments indicate position and orientation. The added faces are returned as a vector. */
std::vector<HMesh::FaceID> create_box(HMesh::Manifold& m, const CGLA::Vec3d& pos, const CGLA::Mat3x3d& _R, double sz);


/** Trace a single face loop starting from h (being a generator of the face loop). touched keeps track of which halfedges have been used. Thus if the
    face loop has already been visited starting from a different generator it will not be */
FaceLoop trace_face_loop(HMesh::Manifold& m, HMesh::HalfEdgeAttributeVector<int>& touched, HMesh::HalfEdgeID h);

/** Creates a vector of all the faceloops of the Manifold m passed as input. */
std::vector<FaceLoop> find_face_loops(HMesh::Manifold& m);

void face_loops_compute_contained_area(HMesh::Manifold& m, std::vector<FaceLoop>& face_loops);


/** Collapses all the face loops in m which are passed as input in the second argument. However, face loops which conain
 the same face twice (i.e. face loops with double cross faces) are not contracted. */
void collapse_face_loops(HMesh::Manifold& m, std::vector<FaceLoop>& face_loops);

/** collapses all double cross quads, i.e. quads that are contained twice in the same face loop. Each collapse splits a face loop in two, and
 the end result is a mesh with no double cross quads. */
int collapse_double_crossed_quads(HMesh::Manifold& m);

/** splits all double cross quads, i.e. quads that are contained twice in the same face loop. Each split of a quad produces two new quads. It works
 by splitting along the diagonal and then inserting a valence two vertex. This also splits the face loop in two, and
 the end result is a mesh with no double cross quads. */
int split_double_crossed_quads(HMesh::Manifold& m);


/** This function takes a set of halfedges (hset) and cuts the mesh along these edges inserting instead a new set of faces along the edges of the cut.
 The new faces are returned as a FaceSet.*/
HMesh::FaceSet extrude_along_edge_loop(HMesh::Manifold& m, const HMesh::HalfEdgeSet& hset);

/** This function extrudes a set of faces. The created face loop is returned as a new face set. Note that the extruded faces are not harmed in the process. */
HMesh::FaceSet extrude_face_set(HMesh::Manifold& m, const HMesh::FaceSet& face_set);

/** Extrude m along a set of halfedges (second argument). The new faces are returned as a FaceSet.  */
HMesh::FaceSet extrude_halfedge_set(HMesh::Manifold& m, HMesh::HalfEdgeSet& halfedge_set);

/** Remove the FacLoop l from m.*/
bool remove_face_loop(HMesh::Manifold& m, const FaceLoop& l);

/** Refine a face loop by splitting it longitudinally.*/
int refine_loop(HMesh::Manifold& m);

/** This function first finds all face loops, sorts them according to the area enclosed on one side of the face loop (the one with the smallest area).
 It then proceeds to remove the face loop with the smallest area. */
void kill_face_loop(HMesh::Manifold& m);

/** Constructs a skeleton from the system of face loops of a mesh m. The face loops are sorted in order of decreasing cylindricity. Then a skeleton
 edge is created for each loop as long as it does not contain a face of a previously seleced face loop. Finally skeleton edges inherit the mesh connectivity
 (roughly speaking) and the skeleton is returned as a graph. */
Geometry::AMGraph3D face_loop_skeleton(HMesh::Manifold& m);



#endif /* face_loop_hmesh */
