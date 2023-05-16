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
    
    /** Class that allows transformations between object coordinates and a voxel grid.
     The object coordinates is simply the coordinate system that contains some data, and
     this class maps them to a coordinate system where voxels are at integer positions.
     While this seems trivial, there are several ways to do it and also pitfalls. This class is
     not meant as a simple solution to the challenge of getting this transformation right but
     as a flexible tool for setting up the transformations. */
    class XForm
    {
        CGLA::Vec3d llf;
        CGLA::Vec3d urt;
        CGLA::Vec3d scale;
        CGLA::Vec3d offset;
        CGLA::Vec3i DIM;
    public:
        XForm() {}
        
        /** Construct an XForm with no margin and where the scale and offset are controlled by the user.
         The _DIM argument is the dimensions of the voxel grid.
         The scale multiplied onto a vector should transform it from the object coordinate scale to the voxel
         grid scale. Likewise, the offset translates a point from object coordinates to voxel grid coordinates.
         */
        XForm(const CGLA::Vec3i& _DIM, double _scale = 1.0, const CGLA::Vec3d& _offset = CGLA::Vec3d(0)):
        llf(0), urt(_DIM-CGLA::Vec3i(1)), scale(_scale), offset(_offset), DIM(_DIM) {}
        
        /** Construct from the corners of the object's bounding volume (lower, left, front) and
            (upper, right, top), the volume dimensions and optionally a margin. */
        XForm(const CGLA::Vec3d& _llf, const CGLA::Vec3d& _urt, const CGLA::Vec3i& _DIM, double margin=0.0):
        llf(_llf), urt(_urt), DIM(_DIM)
        {
            if(urt[0]<llf[0]) {
                scale = CGLA::Vec3d(1.0);
                offset = CGLA::Vec3d(0.0);
                llf = CGLA::Vec3d(0);
                urt = CGLA::Vec3d(DIM-CGLA::Vec3i(1));
            }
            else {
                CGLA::Vec3d d = urt - llf;
                offset = CGLA::Vec3d(d.max_coord() * margin);
                scale = CGLA::Vec3d(DIM-CGLA::Vec3i(1))/(d + 2.0 * offset);
            }
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
        const CGLA::Vec3d inv_scale() const {
            return CGLA::Vec3d(1.0)/scale;
        }
        
        /// Get volume dimensions
        CGLA::Vec3i get_dims() const {
            return DIM;
        }

        /// Return the scale: ratio of voxel size to object size.
        const CGLA::Vec3d& get_scale() const {
            return scale;
        }
        
        /** Return the lower, left, front corner of the bounding box. Note this is in object
        coordinates. */
        const CGLA::Vec3d& get_llf() const {
            return llf;
        }
        
        /** Return the upper, right, top corner of the bounding box. Note this is in object
        coordinates. */
        const CGLA::Vec3d& get_urt() const {
            return urt;
        }

        
    };
}

#endif
