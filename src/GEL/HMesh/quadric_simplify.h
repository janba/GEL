/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file quadric_simplify.h
 * @brief Garland Heckbert simplification.
 */
#ifndef __HMESH_QUADRIC_SIMPLIFY__H
#define __HMESH_QUADRIC_SIMPLIFY__H

#include "Manifold.h"

namespace HMesh
{
    /** \brief Garland Heckbert simplification in our own implementation. 
    keep_fraction is the fraction of vertices to retain. The singular_thresh defines how
    small singular values from the SVD we accept. It is relative to the greatest singular value. 
    If choose_optimal_positions is true, we reposition vertices. Otherwise the vertices are a subset
    of the old vertices. */
    void quadric_simplify(Manifold& m, double keep_fraction, double singular_thresh = 0.0001, bool choose_optimal_positions = true);
    void quadric_simplify(Manifold& m, VertexAttributeVector<int> mask, double keep_fraction, double singular_thresh = 0.0001, bool choose_optimal_positions = true);

}
#endif