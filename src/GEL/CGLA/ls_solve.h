//
//  ls_solve.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 19/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef GEL_ls_solve_h
#define GEL_ls_solve_h

#include <GEL/CGLA/Mat2x2d.h>
#include <GEL/CGLA/Mat2x2f.h>
#include <GEL/CGLA/Mat2x3d.h>
#include <GEL/CGLA/Mat2x3f.h>
#include <GEL/CGLA/Mat3x3d.h>
#include <GEL/CGLA/Mat3x3f.h>
#include <GEL/CGLA/Mat4x4d.h>
#include <GEL/CGLA/Mat4x4f.h>
#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Vec4d.h>
#include <GEL/CGLA/Vec4f.h>

namespace CGLA {
        
    template<class V>
    V ls_solve(const std::vector<V>& A, const std::vector<typename V::ScalarType>& b)
    {
        using M = typename VecT_to_MatT<V>::MatT;
        V ATb(0);
        M ATA(0);
        for(int n = 0; n < A.size(); ++n) {
            for (int i=0; i < ATb.get_dim(); ++i) {
                ATA[i] += A[n] * A[n][i];
            }
            ATb += b[n]*A[n];
        }
        return invert(ATA)*ATb;
    }
    
}

#endif
