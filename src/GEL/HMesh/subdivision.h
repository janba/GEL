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
    void cc_split(Manifold&, Manifold&);
    
    void loop_split(Manifold& m_in, Manifold& m);
    
    void root3_subdivide(Manifold&, Manifold&);
    
    void rootCC_subdivide(Manifold& m_in, Manifold& m);

    void butterfly_subdivide(Manifold& m_in, Manifold& m);

    void cc_smooth(Manifold&);

    void loop_smooth(Manifold&);
}

#endif