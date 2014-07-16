/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file flatten.h
 * @brief Flattening as used in parametrization.
 */

#ifndef __HMESH_FLATTEN_H__
#define __HMESH_FLATTEN_H__

namespace HMesh
{
    class Manifold;

    enum WeightScheme {FLOATER_W, HARMONIC_W, LSCM_W, BARYCENTRIC_W};


    /** \brief This function flattens a mesh with a simple boundary, based on a weight scheme.It is mostly for showing mesh parametrization methods. 
    The current mesh MUST have a SINGLE boundary loop.
    This loop is mapped to the unit circle in a regular fashion (equal angle intervals).
    All non boundary vertices are placed at the origin. Then the system is relaxed iteratively using the weight scheme given as argument. */
    void flatten(Manifold& m, WeightScheme ws);
}

#endif