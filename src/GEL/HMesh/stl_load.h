/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HMesh/stl_load.h
 * @brief Load Manifold from stl file
 */

#ifndef __HMESH_STLLOAD__H__
#define __HMESH_STLLOAD__H__

#include <string>
#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
    /** Load a STL file. 
        The first argument is a string containing the file name (including path) 
     and the second is the Manifold into which the mesh is loaded. The third argument
     is an attribute vector containing the indices of the original
    points.  */
     
    bool stl_load(const std::string&, Manifold& m, VertexAttributeVector<int>& orig_vertex_indices);
    bool stl_load(const std::string&, Manifold& m);
}
#endif
