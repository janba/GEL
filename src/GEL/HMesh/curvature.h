/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file curvature.h
 * @brief Compute various curvature measures from meshes.
 */

#ifndef __MESHEDIT_CURVATURE_H__
#define __MESHEDIT_CURVATURE_H__

#include <vector>
#include "Manifold.h"

namespace CGLA
{
    class Vec3d;
    class Vec2d;
    class Mat2x2d;
    class Mat3x3d;
    class Ma4x4d;
}

namespace HMesh
{
    class Manifold;
    template<typename ITEM>
    class VertexAttributeVector;

    double mixed_area(const Manifold& m,
                        VertexID v);

    double barycentric_area(const Manifold& m, 
                            VertexID v);

    void unnormalized_mean_curvature_normal(const Manifold& m, 
                                            VertexID v, 
                                            CGLA::Vec3d& curv_normal, 
                                            double& w_sum);

    CGLA::Vec3d mean_curvature_normal(  const Manifold& m, 
                                        VertexID v);

    double sum_curvatures(  const Manifold& m, 
                            VertexAttributeVector<double>& curvature);


    double gaussian_curvature_angle_defect( const Manifold& m, 
                                            VertexID v);

    CGLA::Mat3x3d curvature_tensor( const Manifold& m, 
                                    HalfEdgeID h);

    CGLA::Mat3x3d curvature_tensor_from_edge( const Manifold& m, 
                                              HalfEdgeID h);


    void curvature_tensor_paraboloid(   const Manifold& m, 
                                        VertexID v,
                                        CGLA::Mat2x2d& curv_tensor, 
                                        CGLA::Mat3x3d& frame);

    void curvature_tensors_from_edges(  const Manifold& m, 
                                        VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors);

    void smooth_curvature_tensors(  const Manifold& m, 
                                    VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors);

    void gaussian_curvature_angle_defects(  const Manifold& m, 
                                            VertexAttributeVector<double>& curvature, 
                                            int smooth_steps=0);

    void mean_curvatures(   const Manifold& m, 
                            VertexAttributeVector<double>& curvature,
                            int smooth_steps=0);


    void curvature_paraboloids( const Manifold& m, 
                                VertexAttributeVector<CGLA::Vec3d>& min_curv_direction, 
                                VertexAttributeVector<CGLA::Vec3d>& max_curv_direction, 
                               VertexAttributeVector<CGLA::Vec2d>& curvature);


    void curvature_from_tensors(const Manifold& m, 
                                const VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors,
                                VertexAttributeVector<CGLA::Vec3d>& min_curv_direction,
                                VertexAttributeVector<CGLA::Vec3d>& max_curv_direction,
                                VertexAttributeVector<CGLA::Vec2d>& curvature);


}

#endif