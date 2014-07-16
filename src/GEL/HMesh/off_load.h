/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file off_load.h
 * @brief Load Manifold from OFF.
 */

#ifndef __HMESH_OFF_LOAD_HMESH__
#define __HMESH_OFF_LOAD_HMESH__

#include <string>

namespace HMesh
{
    class Manifold;

    /// Load an OFF file (Object File Format). So far, this loader is mostly ensured to load files from the Princeton Shape Benchmark.
    bool off_load(const std::string&, Manifold& m);
}

#endif