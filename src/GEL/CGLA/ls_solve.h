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
    
    template<typename V> class Vec2MatType {};
    template<> struct Vec2MatType<Vec2d> {using Mat = Mat2x2d;};
    template<> struct Vec2MatType<Vec3d> {using Mat = Mat3x3d;};
    template<> struct Vec2MatType<Vec4d> {using Mat = Mat4x4d;};
    template<> struct Vec2MatType<Vec2f> {using Mat = Mat2x2f;};
    template<> struct Vec2MatType<Vec3f> {using Mat = Mat3x3f;};
    template<> struct Vec2MatType<Vec4f> {using Mat = Mat4x4f;};
    
    
    template<class V>
    V ls_solve(const std::vector<V>& A, const std::vector<typename V::ScalarType>& b)
    {
        Vec2MatType<Vec3d>::Mat ATA(0);
        V ATb(0);
        size_t N = A.size();
        for(int n = 0; n < N; ++n){
            ATb += b[n]*A[n];
            for(int i = 0; i < 3; ++i)
                for(int j = 0; j < 3; ++j)
                    ATA[i][j] += A[n][i]*A[n][j];
        }
        return invert(ATA)*ATb;
    }
    
}

#endif
