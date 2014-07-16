/**
 * @file smooth.h
 * @brief Functions for mesh smoothing.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_SMOOTH_H
#define __HMESH_SMOOTH_H

#include <vector>
#include "Manifold.h"

namespace HMesh
{
   enum NormalSmoothMethod {FVM_NORMAL_SMOOTH, BILATERAL_NORMAL_SMOOTH};
    
    /// Simple laplacian smoothing with an optional weight.
    void laplacian_smooth(HMesh::Manifold& m, float t=1.0f, int iter=1);

    /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
    void taubin_smooth(HMesh::Manifold& m, int iter=1);

    /** Smooth meshes by first filtering normals and then refitting the mesh */
    void anisotropic_smooth(HMesh::Manifold& m, int iter, NormalSmoothMethod nsm);

    /// Tangential area weighted smoothing.
    void TAL_smoothing(HMesh::Manifold& m, float w, int iter=1);

}
#endif