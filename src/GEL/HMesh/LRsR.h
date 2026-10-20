/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file LRsR.h
 * @brief Local rotation system reconstruction from a point cloud.
 */

#ifndef HMESH_LRSR_H
#define HMESH_LRSR_H

#include <vector>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/HMesh/Manifold.h>

namespace HMesh::LRsR {

/// Delaunay radius is `radius_factor` times the median nearest-neighbor distance.
/// The 3DDelaunay demo uses 4 for a point cloud.
struct LRsROpts {
    double radius_factor = 4.0;
};

/** Reconstruct a triangle mesh from a point cloud.
 The Delaunay graph, the unoriented adaptive-covariance normals, and the
 tunnel-weighted minimum spanning forest are the same objects the 3DDelaunay
 demo draws. Normals are used only to weight that forest. Each tree is walked
 breadth-first from its lowest-index vertex. At a vertex the link is the
 combinatorial rotation of Delaunay neighbors around that vertex. A closed
 link is a full umbrella and an open link is the partial disk it covers.
 The winding matches triangles chosen at earlier vertices. Where the
 umbrellas put three triangles on one edge, or wind an edge twice the same
 way, the larger triangle is dropped. The triangles are added one face at a
 time, then stitched and cleaned up. */
void reconstruct(const std::vector<CGLA::Vec3d>& points,
                 const LRsROpts& opts,
                 Manifold& result);

}

#endif
