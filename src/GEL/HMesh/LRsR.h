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
 The radius-limited Delaunay graph supplies the candidate triangles. An
 edge under consideration contributes one triangle on the left of its
 normal and one on the right. The candidate minimizes perimeter times area
 divided by the area projected along the edge normal. The edge normal is
 the average of the unoriented adaptive normals at its two vertices. When
 those normals point opposite ways, one is flipped before the average.
 `kGrowthMode` in LRsR.cpp selects the queue. It is 2.
 Mode 1 queues both sides of every tunnel-weighted spanning-tree edge and,
 in the same queue, the free side of each new edge of an accepted triangle.
 Mode 2 starts from the single lowest-score triangle and then queues only
 free sides of edges of accepted triangles. Mode 3 exhausts the tree edges
 first and does not open other edges until that queue is empty. It then
 continues like mode 2 from the triangles already added. Ties break by the
 smaller corner index, then the edge, then the side. The queued score
 multiplies that value by `2 - cos(alpha)` for each face already across one
 of the triangle's edges. Alpha is the angle between the normals once the
 two faces traverse the shared edge in opposite directions. Coplanar faces
 with the same normal keep factor 1. Coplanar faces with opposite normals
 use factor 3. A popped entry whose score has risen is pushed back at the
 current score and is not added yet. Mode 2 sets the queued winding so the
 triangle traverses the proposing edge opposite the face already there, and
 drops a popped triangle that meets any existing neighbor in the same
 direction. A rejected triangle is replaced by
 the next score on the same side of the edge that proposed it.
 A side that already holds a triangle is closed.
 A candidate that would lie on an occupied side of any of its edges is
 skipped for the next score. A candidate that would be the third triangle
 on any edge is skipped the same way. So is a candidate whose wedge at a
 vertex overlaps a triangle already incident there. The edges are projected
 onto the plane perpendicular to the vertex normal. Overlap means a
 non-shared edge of either triangle lies in the other's angle, so a new
 triangle that contains an existing one is skipped along with one that
 falls inside it. An edge shared with the incident triangle is the boundary
 and is allowed. A side with no positive-area Delaunay triangle contributes
 nothing. The two triangles of one edge run in opposite directions along
 that edge. They are both kept only when their wedges do not overlap at
 either endpoint. After every triangle is chosen, each connected component
 is oriented from its first triangle. A neighbor that traverses a shared
 edge in the same direction is reversed, and the walk continues. Only the
 cyclic order of the three vertices is used.
 The return value is how many triangles the queue accepted. Mode 2 does not
 count its first triangle. The run prints `LRsR mode: M` and
 `LRsR iterations: N`. */
int reconstruct(const std::vector<CGLA::Vec3d>& points,
                const LRsROpts& opts,
                Manifold& result);

}

#endif
