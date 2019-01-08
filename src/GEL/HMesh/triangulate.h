/**
 * @file triangulate.h
 * @brief Triangulating the faces of a mesh.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_TRIANGULATE__H
#define __HMESH_TRIANGULATE__H

#include "Manifold.h"

namespace HMesh
{
    enum TriangulationMethod { CLIP_EAR, SHORTEST_EDGE};
    /** Triangulate by connected vertices on the face f.
     The policy indicates if we do ear clip or shortest edge triangulation.
     ear clip is safer, but shortest edge tends to never fail. */
    int triangulate(Manifold& m, FaceID f, TriangulationMethod policy = CLIP_EAR);

    /// Triangulate by connectin
    void triangulate(Manifold& m, TriangulationMethod policy = CLIP_EAR);
}

#endif
