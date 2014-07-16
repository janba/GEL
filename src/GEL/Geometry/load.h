/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Geometry/load.h
 * @brief Load triangle mesh from one of a small number of file formats.
 */

#ifndef __GEOMETRY_LOAD_H__
#define __GEOMETRY_LOAD_H__

#include <string>
#include "TriMesh.h"

namespace Geometry
{
	/// Load a TriMesh from a file. Loader chosen based on extension.
	void load(const std::string &filename, TriMesh &mesh);
}



#endif