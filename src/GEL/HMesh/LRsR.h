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
 The radius-limited Delaunay graph supplies the candidate triangles. Each
 edge of the tunnel-weighted spanning forest contributes one triangle on
 the left of its normal and one on the right. The candidate minimizes
 perimeter times area divided by the area projected along the edge normal.
 The edge normal is the average of the unoriented adaptive
 normals at its two vertices. When those normals point opposite ways, one
 is flipped before the average. A triangle already chosen for an earlier
 forest edge fills the side it lies on when it contains this edge as well,
 and that side is not searched again. A candidate that would lie on an
 occupied side of any of its edges is skipped for the next score. A
 candidate that would be the third triangle on any edge is skipped the
 same way. So is a candidate whose edge at a vertex falls between the
 edges of a triangle already incident there, measured in the plane of the
 vertex normal, unless that edge is shared with the incident triangle.
 A side with no positive-area Delaunay
 triangle contributes nothing. The two new triangles of one edge run in
 opposite directions along that edge. Neighboring edges are not made to agree. */
void reconstruct(const std::vector<CGLA::Vec3d>& points,
                 const LRsROpts& opts,
                 Manifold& result);

}

#endif
