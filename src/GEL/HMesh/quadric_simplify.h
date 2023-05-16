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

#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
    /** \brief Garland Heckbert simplification in our own implementation. 
    keep_fraction is the fraction of vertices to retain. The singular_thresh controls sensitivity to subtle sharp edges and corners. If the
    parameter is close to 0 subtler features are preserved. The err_thresh is a threshold on the quadric error measure itself. The mesh
    will be simplified until keep_fraction is reached, unless the error exceeds err_thresh before that happens. */
    void quadric_simplify(Manifold& m, double keep_fraction, double singular_thresh = 0.0001, double err_thresh=0.0);
}
#endif
