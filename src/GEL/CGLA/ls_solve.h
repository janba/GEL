//
//  ls_solve.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 19/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef GEL_ls_solve_h
#define GEL_ls_solve_h

#include "Mat2x2d.h"
#include "Mat2x2f.h"
#include "Mat2x3d.h"
#include "Mat2x3f.h"
#include "Mat3x3d.h"
#include "Mat3x3f.h"
#include "Mat4x4d.h"
#include "Mat4x4f.h"
#include "Vec2d.h"
#include "Vec2f.h"
#include "Vec3d.h"
#include "Vec3f.h"
#include "Vec4d.h"
#include "Vec4f.h"

namespace CGLA {
        
    template<class V>
    V ls_solve(const std::vector<V>& A, const std::vector<typename V::ScalarType>& b)
    {
        typename VecT_to_MatT<V>::MatT ATA(0);
        V ATb(0);
        for(int n = 0; n < A.size(); ++n) {
            ATA += outer_product(A[n], A[n]);
            ATb += b[n]*A[n];
        }
        return invert(ATA)*ATb;
    }
    
}

#endif
