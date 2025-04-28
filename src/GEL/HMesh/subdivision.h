/**
 * @file subdivision.h
 * @brief Functions for mesh subdivision. Catmull Clark to be precise.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_SUBDIVIDE_H__
#define __HMESH_SUBDIVIDE_H__

namespace HMesh
{
    class Manifold;
    /** Perform a Catmull-Clark split, i.e. a split where each face is divided
    into new quadrilateral faces formed by connecting a corner with a
    point on each incident edge and a point at the centre of the face. */
    void cc_split(Manifold&);
    

    /** Perform a loop style split. The input manifold is assumed to be a triangle
     mesh and each face is split into four faces. */
    void loop_split(Manifold&);
    
    /** Perform one step of Kobbelt's sqrt-3 subdivision. */
    void root3_subdivide(Manifold&);
    
    
    void rootCC_subdivide(Manifold& m);

    /** Perform one step of butterfly subdivision, i.e. interpolatory subdivision
     on triangle meshes. */
    void butterfly_subdivide(Manifold& m);

    /** Perform Catmull Clark smoothing. If this function follows cc_split it results
     in Catmull-Clark subdivison implemented in a factored fashion. */
    void cc_smooth(Manifold&);

    /** This function is inspired by Taubin smoothing and works roughly like CC smooth
     but causes less shrinkage. Very homegrown - use at your own risk. */
    void volume_preserving_cc_smooth(Manifold& m, int iter=1);

    /** Perform Loop smoothing. If this function follows loop_split, the result is
     Loop subdivision implemented in a factored fashion. */
    void loop_smooth(Manifold&);
}

#endif
