/**
 * @file refine_edges.h
 * @brief Tools for subdividing edges of a Manifold.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_REFINE_EDGES_H__
#define __HMESH_REFINE_EDGES_H__

namespace HMesh
{
    class Manifold;

    /// Return the average edge length
    float average_edge_length(const Manifold& m);

    /** Split all edges in mesh passed as first argument which are longer
    than the threshold (second arg) length. A split edge
    results in a new vertex of valence two.*/
    int refine_edges(Manifold& m, float t);
}

#endif