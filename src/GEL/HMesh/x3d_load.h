/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file x3d_load.h
 * @brief Load from X3D
 */

#ifndef __HMESH_X3DLOAD_H__
#define __HMESH_X3DLOAD_H__

#include <string>

namespace HMesh
{
    class Manifold;
    /// Load a mesh from an X3D file. It handles arbitrary polygons.
    bool x3d_load(const std::string& filename, Manifold& m);
}
#endif