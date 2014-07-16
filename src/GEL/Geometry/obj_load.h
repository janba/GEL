/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Geometry/obj_load.h
 * @brief Load a Wavefront OBJ file into a triangle mesh.
 */

#ifndef __GEOMETRY_TRIMESH_FUNCTIONS_H__
#define __GEOMETRY_TRIMESH_FUNCTIONS_H__

#include <string>
#include <vector>
#include "TriMesh.h"

namespace Geometry
{
	/// Load a TriMesh from an OBJ file
	void obj_load(const std::string &filename, TriMesh &mesh);

  /// Load materials from an MRL file
  void mtl_load(const std::string& filename, std::vector<Material>& materials);
}

#endif
