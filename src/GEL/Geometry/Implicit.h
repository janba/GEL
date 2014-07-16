//
//  Implicit.h
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __PointReconstruction__Implicit__
#define __PointReconstruction__Implicit__

#include "../CGLA/Vec3d.h"
#include "RGrid.h"
#include "XForm.h"

namespace Geometry
{
    class Implicit
    {
    public:
        virtual ~Implicit() {}
        virtual double eval(const CGLA::Vec3d& p) const = 0;
        virtual CGLA::Vec3d grad(const CGLA::Vec3d& p) const = 0;
        void push_to_surface(CGLA::Vec3d& p, double tau=0, double max_dist=FLT_MAX) const;
    };
    
    XForm grid_sample(const Implicit& imp, const CGLA::Vec3d& llf, const CGLA::Vec3d& urt,
                      Geometry::RGrid<float>& grid);
    
    class VolumetricImplicit: public Implicit {
        const XForm xform;
        const Geometry::RGrid<float>& grid;
    public:
        VolumetricImplicit(const XForm& _xform, const Geometry::RGrid<float>& _grid):
        xform(_xform), grid(_grid) {}
        virtual double eval(const CGLA::Vec3d& p) const;
        virtual CGLA::Vec3d grad(const CGLA::Vec3d& p) const;
    };
}
#endif /* defined(__PointReconstruction__Implicit__) */
