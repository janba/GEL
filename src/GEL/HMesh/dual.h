/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file dual.h
 * @brief Compute the dual of a mesh.
 */

#ifndef __HMESH_DUAL_H__
#define __HMESH_DUAL_H__

namespace HMesh
{
    class Manifold;

    void dual(Manifold& m);
}

#endif