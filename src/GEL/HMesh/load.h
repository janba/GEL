/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HMesh/load.h
 * @brief Load a Manifold from various types of files.
 */



#ifndef __HMESH_LOAD__H__
#define __HMESH_LOAD__H__

#include <string>

namespace HMesh
{
    class Manifold;
    
    /// Load a geometry file. This could be a PLY, OBJ, X3D, or OFF file
    bool load(const std::string&, Manifold&, bool safe=true);
    
    /** \brief Safely build a manifold.
     The arguments are the number of vertices (no_vertices),  the vector of vertices (vertvec),
     the number of faces (no_faces), a pointer to an array of double values (vert_vec) and an array
     of indices (indices).
     Note that each vertex is three double precision floating point numbers.
     The indices vector is one long list of all vertex indices. Note also that this function
     does not assume that the mesh is manifold. Each face is created independently and then stitched
     along the boundary with adjacent faces. This stitching will fail in non-manifold situations, but
     loading should always complete.  */
    void safe_build(Manifold& m, size_t no_vertices,
                    const double* vertvec,
                    size_t no_faces,
                    const int* facevec,
                    const int* indices);
}

#endif
