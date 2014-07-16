/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file off_save.h
 * @brief Save Manifold to OFF.
 */

#ifndef __HMESH_OFF_SAVE_H__
#define __HMESH_OFF_SAVE_H__

#include <string>

namespace HMesh
{
    class Manifold;
    /** Save in OFF format. */
    bool off_save(const std::string&, HMesh::Manifold& m);

}

#endif
