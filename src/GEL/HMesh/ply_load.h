/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HMesh/ply_load.h
 * @brief Load Manifold from PLY format.
 */

#ifndef __HMESH_PLYLOAD__H__
#define __HMESH_PLYLOAD__H__

#include <string>

namespace HMesh
{
    class Manifold;
    /** Load an Wavefront OBJ file. This is just a simple frontend for the 
    Geometry::obj_load function which loads OBJ files into triangle meshes. 
    Consequently, quads are unfortunately converted to triangles when this
    loader is used. */
    bool ply_load(const std::string&, Manifold& m);
}
#endif