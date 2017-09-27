//
//  XForm.h
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#ifndef PointReconstruction_XForm_h
#define PointReconstruction_XForm_h

#include "../CGLA/Vec3d.h"

namespace Geometry
{
    
    /** Class that allows transformations between a voxel space and object coordinates */
    class XForm
    {
        CGLA::Vec3d llf;
        CGLA::Vec3d urt;
        double scale, margin;
        CGLA::Vec3i DIM;
    public:
        XForm() {}
        
        /** Construct from the corners of the object's bounding volume (lower, left, front) and
            (upper, right, top), the volume dimensions and optionally margin. */
        XForm(const CGLA::Vec3d& _llf, const CGLA::Vec3d& _urt, const CGLA::Vec3i& _DIM, double _margin=0.1):
        llf(_llf), urt(_urt), DIM(_DIM)
        {
            if(urt[0]<llf[0])
            {
                margin = 0.0;
                scale = 1;
                llf = CGLA::Vec3d(0);
                urt = CGLA::Vec3d(DIM);
            }
            else
            {
                CGLA::Vec3d d = urt - llf;
                margin = d.max_coord() * _margin;
                double scale_x = (DIM[0])/(d[0] + 2.0 * margin);
                double scale_y = (DIM[1])/(d[1] + 2.0 * margin);
                double scale_z = (DIM[2])/(d[2] + 2.0 * margin);
                scale = fmin(fmin(scale_x, scale_y), scale_z);
            }
        }
        
        /// Get volume dimensions
        CGLA::Vec3i get_dims() const {return DIM;}
        
        /// Apply: transform from object to voxel coords.
        const CGLA::Vec3d apply(const CGLA::Vec3d& p) const
        {
            return scale*(p-llf+CGLA::Vec3d(margin));
        }
        
        /// Apply inverse: transform from voxel to object coords.
        const CGLA::Vec3d inverse(const CGLA::Vec3d& p) const
        {
            return p/scale + llf - CGLA::Vec3d(margin);
        }
        
        /// Return the inverse scale: ratio of object to voxel size.
        double inv_scale() const {return 1.0/scale;}

        /// Return the scale: ratio of voxel size to object size.
        double get_scale() const {return scale;}
        
        
        /// Printout information about the transformation
        void print()
        {
            std::cout << scale << " " << (-llf+CGLA::Vec3d(margin)) << std::endl;
        }
    };
}

#endif
