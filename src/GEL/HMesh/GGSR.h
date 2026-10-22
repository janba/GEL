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
/// `max_edges` is the Dijkstra hop limit. The default is 10.
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
 `GGSR mst edges: M`. Face finding still walks the full Gabriel graph.
 From each vertex, unit-weight Dijkstra grows a tree of paths at most
 `max_edges` edges long. Each other Gabriel edge between two reached
 vertices closes a face with the tree path that joins its ends. A face is
 stored under the hash of its vertex indices in increasing order, and a
 later discovery of the same vertex set is ignored. The stored boundary
 and starting vertex belong to the first discovery. Two faces that share
 more than three edges drop the larger one. Two that share exactly two
 edges project the vector from each starting vertex to that face's
 barycenter into the plane perpendicular to the normal at that vertex. A
 positive dot product of the projected vectors drops the larger face as
 well. The larger face has more edges, or, when the counts are equal, the
 greater Euclidean perimeter. A pair of equal size is kept, as is every
 pair that shares one or three edges. The kept polygons are added to
 `result`, stitched, and each connected set of faces is oriented so a
 shared edge is traversed in opposite directions. The return value is the
 number of kept faces. It is also printed as `GGSR faces: N`. */
int reconstruct(const std::vector<CGLA::Vec3d>& points,
                const GGSROpts& opts,
                Manifold& result);

}

#endif
