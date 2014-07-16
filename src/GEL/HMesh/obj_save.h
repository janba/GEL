/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file obj_save.h
 * @brief Save Manifold to OBJ.
 */

#ifndef __HMESH_OBJSAVE__H__
#define __HMESH_OBJSAVE__H__

#include <string>

namespace HMesh
{
    class Manifold;
    /// \brief Save in Wavefront OBJ format. 
    bool obj_save(const std::string&, Manifold& m);

}
#endif