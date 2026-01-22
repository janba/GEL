/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file stl_save.h
 * @brief Save Manifold to STL.
 */

#ifndef __HMESH_STLSAVE__H__
#define __HMESH_STLSAVE__H__

#include <string>
#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
    class Manifold;
    /// \brief Save in STL format. 
    bool stl_save(const std::string&, Manifold& m, bool is_binary = false);

}
#endif
