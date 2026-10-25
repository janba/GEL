/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file JohnsonCycles.h
 * @brief Length-bounded elementary cycles of an undirected graph.
 */

#ifndef GEOMETRY_JOHNSON_CYCLES_H
#define GEOMETRY_JOHNSON_CYCLES_H

#include <vector>
#include <GEL/Geometry/Graph.h>

namespace Geometry {

/// One elementary cycle. `vertices` walks the cycle once, starting at the
/// least vertex id and leaving it toward the smaller of its two neighbors on
/// the cycle. `start` is that least vertex.
struct ElementaryCycle {
    std::vector<AMGraph::NodeID> vertices;
    AMGraph::NodeID start = AMGraph::InvalidNodeID;
};

/**
 * Elementary cycles of the undirected graph `graph`.
 *
 * An elementary cycle is a closed walk of at least three edges that repeats
 * no vertex. Each such cycle with at most `max_edges` edges is returned
 * once. The two directions of a cycle are the same cycle.
 *
 * The enumeration is Johnson's circuit search (SIAM J. Comput. 4, 1975),
 * length-bounded and run on an undirected graph. Vertices are processed in
 * increasing id. The search from a vertex keeps only vertices of id at least
 * that large, so every cycle is grown from its least vertex. A vertex entered
 * on a path of some length stays closed to longer paths until a cycle shows
 * that it can still return to the start inside the length bound.
 *
 * `max_edges < 3` yields no cycles.
 */
std::vector<ElementaryCycle> elementary_cycles(const AMGraph& graph, int max_edges);

}

#endif
