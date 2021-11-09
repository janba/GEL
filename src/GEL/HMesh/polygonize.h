//
//  polygonize.h
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __PointReconstruction__polygonize__
#define __PointReconstruction__polygonize__

#include <iostream>

#include <GEL/Geometry/XForm.h>
#include <GEL/Geometry/RGrid.h>
#include <GEL/HMesh/Manifold.h>

namespace HMesh
{
    float clamp_interpolate(const Geometry::RGrid<float>& grid, const CGLA::Vec3d& v);
    CGLA::Vec3f clamp_trilin_grad(const Geometry::RGrid<float>& grid, const CGLA::Vec3d& v);

    /**
     @brief Computes a polygonal mesh from a volumetric isocontour
     @param xform is the transformation from voxel indices back to the domain of the implicit being polygonized.
     @param grid is the voxel grid.
     @param mani is the manifold into which the output mesh is written.
     @param tau is the threshold (or isovalue)
     @param make_triangles tells whether we output triangles (default: true) or quads (false)
     @param high_is_inside tells whether values greater than @tau (defalut: true) are interior or exterior
     @param pre_smooth_steps is how many steps of (Taubin) smoothing to apply before projecting vertices onto surface.
     @returns Nothing
     
     This function computes an iso surface using the method of dual contouring. For each voxel that is inside,
     it visits the six neighbors and outputs a quad if that neighbor is outside. This leads to a cuberille mesh.
     Afterwards, vertices are smoothed and finally placed on the isocontour by using a line search in the
     direction of the vertex normal. Dual contouring is very simple and leads to bette triangles than marching
     cubes. On the flip side, the vertex placement is arguably a bit more ad hoc.
     */
    void volume_polygonize(const Geometry::XForm& xform, const Geometry::RGrid<float>& grid,
                           HMesh::Manifold& mani, float tau, bool make_triangles=true, bool high_is_inside=true,  int pre_smooth_steps=0 );
}

#endif /* defined(__PointReconstruction__polygonize__) */

