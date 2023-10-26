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
#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
   enum NormalSmoothMethod {FVM_NORMAL_SMOOTH, BILATERAL_NORMAL_SMOOTH};
    
    inline CGLA::Vec3d laplacian(const Manifold& m, VertexID v)
    {
        CGLA::Vec3d p(0);
        int cnt = 0;
        for(auto vn: m.incident_vertices(v)) {
            p += m.pos(vn);
            ++cnt;
        }
        return p / cnt - m.pos(v);
    }

    CGLA::Vec3d cot_laplacian(const Manifold& m, VertexID v);
    
    /// Simple laplacian smoothing with an optional weight.
    void laplacian_smooth(HMesh::Manifold& m, float t=1.0f, int iter=1);

    /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
    void taubin_smooth(HMesh::Manifold& m, int iter=1);

    /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage. This version uses the cotan weights.
    void taubin_smooth_cot(HMesh::Manifold& m, int iter=1);

    /** Smooth meshes by first filtering normals and then refitting the mesh */
    void anisotropic_smooth(HMesh::Manifold& m, int iter, NormalSmoothMethod nsm);

    /// Tangential area weighted smoothing.
    void TAL_smoothing(HMesh::Manifold& m, float w, int iter=1);

    /// Smooth mesh by regularizing quads, i.e. making them more rectangular.
    void regularize_quads(HMesh::Manifold &m, float w, float shrink, int max_iter);
}
#endif
