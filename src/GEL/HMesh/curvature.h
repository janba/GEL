/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file curvature.h
 * @brief Compute various curvature measures from meshes.
 */

#ifndef MESHEDIT_CURVATURE_H
#define MESHEDIT_CURVATURE_H

#include <vector>
#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
    class Manifold;
    template<typename ITEM>
    class VertexAttributeVector;

    /**
     * @brief Smooth a vector field on the mesh vertices.
     * @param m The input mesh.
     * @param vec The vector field to smooth (in-place).
     * @param smooth_steps Number of smoothing iterations.
     */
    void smooth_vectors_on_mesh(const Manifold& m, VertexAttributeVector<CGLA::Vec3d>& vec, int smooth_steps);

    /**
     * @brief Compute the mixed (approximate Voronoi) area associated with a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The mixed area for the vertex.
     */
    double mixed_area(const Manifold& m, VertexID v);

    /**
     * @brief Compute the barycentric area (1/3 of the incident triangles' areas) associated with a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The barycentric area for the vertex.
     */
    double barycentric_area(const Manifold& m, 
                            VertexID v);

    /**
     * @brief Compute the mean curvature normal at a vertex using the cotan formula.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The mean curvature normal vector.
     */
    CGLA::Vec3d mean_curvature_normal(  const Manifold& m, 
                                        VertexID v);

    /**
     * @brief Compute the mean curvature at a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The mean curvature value computed as the ratio of the length of the mean curvature normal to the mixed area.
     */
    double mean_curvature(const Manifold& m, VertexID v);

    /**
     * @brief Compute the angle defect at a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The angle defect (2pi - sum of incident angles).
     */
    double angle_defect(const Manifold& m, VertexID v);

    /**
     * @brief Compute the Gaussian curvature at a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return The Gaussian curvature value computed as the ratio of the angle defect to the mixed area.
     */
    double gaussian_curvature(const Manifold& m, VertexID v);

        /**
     * @brief Compute Gaussian curvature for all vertices.
     * @param m The input mesh.
     * @param curvature Output: Gaussian curvature per vertex.
     * @param smooth_steps Number of smoothing iterations (default 0).
     */
    void gaussian_curvature(const Manifold& m, 
        VertexAttributeVector<double>& curvature, 
        int smooth_steps=0);

    /**
     * @brief Compute mean curvature for all vertices.
     * @param m The input mesh.
     * @param curvature Output: mean curvature per vertex.
     * @param smooth_steps Number of smoothing iterations (default 0).
     */
    void mean_curvature(const Manifold& m, 
        VertexAttributeVector<double>& curvature,
        int smooth_steps=0);

    /**
     * @brief Struct holding principal curvature values and directions at a vertex.
     */
    struct PrincipalCurvatures {
        CGLA::Vec3d min_curv_direction;  ///< Direction of minimum curvature
        CGLA::Vec3d max_curv_direction;  ///< Direction of maximum curvature
        double min_curvature;        ///< Minimum curvature
        double max_curvature;        ///< Maximum curvature
    };

    /**
     * @brief Compute the principal curvatures and directions at a vertex.
     * @param m The input mesh.
     * @param v The vertex ID.
     * @return PrincipalCurvatures struct with directions and values.
     * 
     * Note that this function computes the principal curvatures by locally fitting
     * a paraboloid to the 1-ring vertices of the mesh. Based on this a 2x2 matrix
     * known as the shape operator is computed.
     */
    PrincipalCurvatures principal_curvatures( const Manifold& m, VertexID v);


    /**
     * @brief Compute curvature tensors for all vertices using edge-based method.
     * @param m The input mesh.
     * @param curvature_tensors Output: curvature tensor per vertex.
     * 
     * This function computes curvature tensors directly from edges using the bending
     * of the mesh to compute a 3x3 matrix that will have an eigenvector in the edge 
     * direction - perpendicular to the direction of curvature. If we construct this
     * matrix for all incident edges of a vertex, we obtain a matrix whose eigenvalues
     * and eigenvectors correspond to principal curvatures and curvature directions 
     * although the index of the minimum direction corresponds to the maximum curvature 
     * and vice versa.
     */
    void curvature_tensors_from_edges(const Manifold& m, 
        VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors);

    /**
     * @brief Smooth curvature tensors on the mesh.
     * @param m The input mesh.
     * @param curvature_tensors Curvature tensors to smooth (in-place).
     */
    void smooth_curvature_tensors(const Manifold& m, 
        VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors);

    /**
     * @brief Extract principal curvatures and directions from curvature tensors.
     * @param m The input mesh.
     * @param curvature_tensors Input: curvature tensor per vertex.
     * @param min_curv_direction Output: direction of minimum curvature per vertex.
     * @param max_curv_direction Output: direction of maximum curvature per vertex.
     * @param curvature Output: (min, max) curvature values per vertex.
     */
    void curvature_from_tensors(const Manifold& m, 
        const VertexAttributeVector<CGLA::Mat3x3d>& curvature_tensors,
        VertexAttributeVector<CGLA::Vec3d>& min_curv_direction,
        VertexAttributeVector<CGLA::Vec3d>& max_curv_direction,
        VertexAttributeVector<CGLA::Vec2d>& curvature);

}

#endif
