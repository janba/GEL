/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HMesh/obj_load.h
 * @brief Load Manifold from OBJ file
 */

#ifndef __HMESH_OBJLOAD__H__
#define __HMESH_OBJLOAD__H__

#include <string>
#include "Manifold.h"

namespace HMesh
{
    /** Load a Wavefront OBJ file. 
        The first argument is a string containing the file name (including path) 
     and the second is the Manifold into which the mesh is loaded. The third and
     final argument is a boolean which indicates whether safe loading is used. This 
     argument defaults to true. Safe loading means that all the faces are loaded
     individually and then stitched in a second pass. 
     If safe is false, the faces are loaded and stitched at the same time, and this
     procedure cannot handle non-manifold situations. */
     
    bool obj_load(const std::string&, Manifold& m, VertexAttributeVector<int>& orig_vertex_indices);
    bool obj_load(const std::string&, Manifold& m);
}
#endif
