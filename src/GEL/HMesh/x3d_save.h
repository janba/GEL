/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file x3d_save.h
 * @brief Save to X3D
 */
#ifndef __HMESH_X3D_SAVE_H__
#define __HMESH_X3D_SAVE_H__

#include <string>

namespace HMesh
{
    class Manifold;
    /// Save mesh to x3d file.
    bool x3d_save(const std::string&, Manifold& m);
}
#endif