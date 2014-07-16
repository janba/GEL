/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Geometry/ply_load.h
 * @brief Load files of the PLY format.
 */

#ifndef __GEOMETRY_PLYLOAD_H_
#define __GEOMETRY_PLYLOAD_H_


#include "../Geometry/TriMesh.h"

namespace Geometry
{
	/** Load geometry from a ply file into a TriMesh. This is a very crude loader which only extracts the 
		raw geometry. */
	bool ply_load(const std::string& fn, Geometry::TriMesh& mesh);
}

#endif