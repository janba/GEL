/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file LocalDijkstraFaces.h
 * @brief Faces formed by chords of a hop-limited Dijkstra tree.
 */

#ifndef GEOMETRY_LOCAL_DIJKSTRA_FACES_H
#define GEOMETRY_LOCAL_DIJKSTRA_FACES_H

#include <vector>
#include <GEL/Geometry/Graph.h>

namespace Geometry {

/// One face closed by a chord of a hop-limited Dijkstra tree.
/// `vertices` is the boundary in the order of that first discovery.
/// `start` is the vertex the search grew from.
struct DijkstraFace {
    std::vector<AMGraph::NodeID> vertices;
    AMGraph::NodeID start = AMGraph::InvalidNodeID;
};

/**
 * Faces from hop-limited Dijkstra trees on `graph`.
 *
 * Each edge has weight 1. From every vertex, Dijkstra therefore grows a
 * breadth-first tree and stops when a vertex is `max_edges` edges from the
 * root. Vertices farther away are not reached. A graph edge that is not in
 * that tree, and whose ends were both reached, closes one face: the edge
 * plus the unique tree path between its ends.
 *
 * Faces are kept in a hash table. The key is the vertex indices sorted into
 * increasing order, and the hash is computed from that sequence. A later
 * discovery of the same vertex set is ignored, so the stored boundary and
 * starting vertex belong to the first discovery.
 *
 * `max_edges < 1` yields no faces. A cycle has at least three vertices.
 */
std::vector<DijkstraFace> local_dijkstra_faces(const AMGraph3D& graph, int max_edges);

}

#endif
