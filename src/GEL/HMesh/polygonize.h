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

#include "../Geometry/XForm.h"
#include "../Geometry/RGrid.h"
#include "Manifold.h"

namespace HMesh
{
    void polygonize(const Geometry::XForm& xform,
                    const Geometry::RGrid<float>& grid,
                    std::vector<CGLA::Vec3d>& quad_vertices, float tau=0);
    
    void volume_polygonize(const Geometry::XForm& xform, const Geometry::RGrid<float>& grid,
                           HMesh::Manifold& mani, float tau);
}

#endif /* defined(__PointReconstruction__polygonize__) */

