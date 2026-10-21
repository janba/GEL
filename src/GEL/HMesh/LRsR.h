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
 Triangles are added smallest score first. The queue starts with both
 sides of every tunnel-weighted forest edge. Adding a triangle brings each
 of its other edges into consideration when that edge is not already
 there, and queues the other side of the new edge. A rejected triangle is
 replaced by the next score on the same side of the edge that proposed it.
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
 The return value is how many accepted triangles were proposed by an edge
 outside the forest. The count is also printed as `LRsR iterations: N`. */
int reconstruct(const std::vector<CGLA::Vec3d>& points,
                const LRsROpts& opts,
                Manifold& result);

}

#endif
