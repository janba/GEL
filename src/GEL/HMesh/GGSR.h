/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file GGSR.h
 * @brief Gabriel graph surface reconstruction.
 */

#ifndef HMESH_GGSR_H
#define HMESH_GGSR_H

#include <vector>
#include <GEL/CGLA/Vec.h>
#include <GEL/HMesh/Manifold.h>

namespace HMesh::GGSR {

/// Gabriel search radius is `radius_factor` times the median nearest-neighbor
/// distance. That absolute radius is passed to `Geometry::gabriel_edges`.
/// `max_edges` is the maximum number of edges on a traced face. The default is 10.
struct GGSROpts {
    double radius_factor = 4.0;
    int max_edges = 10;
};

/** Find Gabriel-graph faces and write them into `result`.
 The point cloud is connected by the radius-limited Gabriel graph. Vertex
 normals are the unoriented adaptive covariance estimate. The
 tunnel-weighted minimum spanning forest of that graph is then built, one
 tree per connected component. An edge of length `l` costs
 `l * (1 + |ê·n| + |ê·m|)`. The number of forest edges is printed as
 `GGSR mst edges: M`. Face finding walks the full Gabriel graph.
 At each vertex, every pair of incident edges is a corner. Walks are
 seeded from those corners and traced from both edges. The cross product
 of the two edges sets the normal, with the corner to the left of travel.
 At each later vertex every other incident edge is a candidate. A positive
 turn back to an earlier vertex of the same walk closes that shorter
 cycle. Otherwise a positive turn whose far vertex is adjacent to the
 start is taken, since the cycle then closes on the next step. The
 shortest of those cycles wins. Otherwise the greatest signed turn in the
 plane perpendicular to the running normal is taken. Equal turns keep the
 smaller neighbor id. The normal is then replaced by the accepted corner's
 cross product, flipped so it stays in the same hemisphere. Tracing stops
 when the walk returns to its start or the next edge would pass
 `max_edges`. A finished cycle is simple and has at least three edges.
 The tracer can claim each corner of that cycle and refuse to start from
 or step onto a claimed corner. That single-use check is currently
 switched off (`kExclusiveCorners` in `GGSR.cpp`), so one corner may
 belong to several cycles. The same undirected cycle is kept once,
 starting at its least vertex id and leaving toward the smaller cycle
 neighbor. Not every cycle of the graph is found. The polygons are added
 to `result`, stitched, and each connected set of faces is oriented so a
 shared edge is traversed in opposite directions. The return value is the
 number of cycles. It is also printed as `GGSR faces: N`. */
int reconstruct(const std::vector<CGLA::Vec3d>& points,
                const GGSROpts& opts,
                Manifold& result);

}

#endif
