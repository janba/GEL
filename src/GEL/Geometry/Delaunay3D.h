/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Delaunay3D.h
 * @brief Radius-limited 3D Delaunay edges.
 */

#ifndef GEOMETRY_DELAUNAY3D_H
#define GEOMETRY_DELAUNAY3D_H

#include <vector>
#include <GEL/Geometry/Graph.h>

namespace Geometry {

/**
 * Delaunay edges of length at most `radius`.
 *
 * The cell around each site starts as a cube of half-extent slightly
 * larger than `radius / 2`, centered on the site. That cube is the
 * bootstrap: its eight corners are the first Voronoi vertices, and each
 * corner stores the three cube planes that meet there. The six cube
 * planes are artificial. They never become graph edges.
 *
 * Neighbors are the other sites returned by the kD-tree at distance at
 * most `radius`, clipped nearest first. Each neighbor adds the
 * perpendicular bisector through the midpoint. A vertex outside that
 * plane is cut off. Each edge into such a vertex is the intersection of
 * two earlier planes, so the new vertex is the intersection of those two
 * planes with the bisector.
 *
 * A neighbor is kept only when its bisector still meets the finished
 * cell in a polygon of positive area. An edge longer than `radius` is
 * outside the search. A Delaunay edge of length at most `radius` whose
 * Voronoi face misses the cube is not certified by this neighborhood and
 * is omitted: the empty ball can sit farther than `radius / 2` from the
 * site, which happens for slivers.
 *
 * Cells are independent and are clipped on C++ threads. The returned
 * graph has one node per input point, in the same order. `radius <= 0`
 * yields those nodes and no edges.
 */
AMGraph3D delaunay_edges(const std::vector<CGLA::Vec3d>& points, double radius);

/**
 * Gabriel edges of length at most `radius`.
 *
 * An edge is Gabriel when the closed ball with that edge as diameter
 * contains no other input point. Every such edge has length at most
 * `radius`, and any point that could lie in its diametral ball is within
 * `radius` of both ends, so the same kD-tree neighborhood decides it.
 * The test is the bisector of each other neighbor at the edge midpoint:
 * the midpoint has to lie strictly inside every one of those halfspaces.
 * No Voronoi cell is built.
 *
 * The search is parallel over sites. The graph has one node per input
 * point, in the same order. `radius <= 0` yields those nodes and no edges.
 */
AMGraph3D gabriel_edges(const std::vector<CGLA::Vec3d>& points, double radius);

}

#endif
