//
//  XForm.h
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#ifndef PointReconstruction_XForm_h
#define PointReconstruction_XForm_h

#include <GEL/CGLA/Vec3d.h>

namespace Geometry
{
    
    /** Class that allows transformations between object coordinates and a voxel grid. */
    class XForm
    {
        CGLA::Vec3d llf;
        CGLA::Vec3d urt;
        double scale;
        CGLA::Vec3d offset;
        CGLA::Vec3i DIM;
    public:
        XForm() {}
        
        /** Construct an XForm suitable for when the voxels (i.e. grid points) are at half integeger coordinates.*/
        XForm(const CGLA::Vec3i& _DIM, double _scale = 1.0, const CGLA::Vec3d&  _offset = CGLA::Vec3d(-0.5)):
        llf(0), urt(_DIM), scale(_scale), offset(_offset), DIM(_DIM) {}
        
        /** Construct from the corners of the object's bounding volume (lower, left, front) and
            (upper, right, top), the volume dimensions and optionally margin. */
        XForm(const CGLA::Vec3d& _llf, const CGLA::Vec3d& _urt, const CGLA::Vec3i& _DIM, double _margin=0.1):
        llf(_llf), urt(_urt), DIM(_DIM)
        {
            if(urt[0]<llf[0]) {
                scale = 1;
                offset = CGLA::Vec3d(0.0);
                llf = CGLA::Vec3d(0);
                urt = CGLA::Vec3d(DIM-CGLA::Vec3i(1));
            }
            else {
                CGLA::Vec3d d = urt - llf;
                double margin = d.max_coord() * _margin;
                double scale_x = (DIM[0]-1)/(d[0] + 2.0 * margin);
                double scale_y = (DIM[1]-1)/(d[1] + 2.0 * margin);
                double scale_z = (DIM[2]-1)/(d[2] + 2.0 * margin);
                scale = fmin(fmin(scale_x, scale_y), scale_z);
                offset = CGLA::Vec3d(margin);
            }
        }
        
        /// Get volume dimensions
        CGLA::Vec3i get_dims() const {
            return DIM;
        }
        
        /// Apply: transform from object to voxel coords.
        const CGLA::Vec3d apply(const CGLA::Vec3d& p) const {
            return scale*(p - llf + offset);
        }
        
        /// Apply inverse: transform from voxel to object coords.
        const CGLA::Vec3d inverse(const CGLA::Vec3d& p) const {
            return p/scale + llf - offset;
        }
        
        /// Return the inverse scale: ratio of object to voxel size.
        double inv_scale() const {
            return 1.0/scale;
        }

        /// Return the scale: ratio of voxel size to object size.
        double get_scale() const {
            return scale;
        }
        
    };
}

#endif
