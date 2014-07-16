/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HMesh/load.h
 * @brief Load a Manifold from various types of files.
 */



#ifndef __HMESH_LOAD__H__
#define __HMESH_LOAD__H__

#include <string>

namespace HMesh
{
    class Manifold;

    /// Load a geometry file. This could be a PLY, OBJ, X3D, or OFF file
    bool load(const std::string&, Manifold& m);
}

#endif