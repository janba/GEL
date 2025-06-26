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
    V lin_solve(const typename VecT_to_MatT<V>::MatT& Ap, const V& b)
    {
        using M = typename VecT_to_MatT<V>::MatT;
        using S = typename V::ScalarType;
        M A = Ap;
        V x = b;
        size_t dim = A.get_v_dim();

        // Forward elimination with partial pivoting
        for (size_t i = 0; i < dim; ++i) {
            // Pivoting
            size_t max_row = i;
            for (size_t k = i + 1; k < dim; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                    max_row = k;
            }
            std::swap(A[i], A[max_row]);
            std::swap(x[i], x[max_row]);

            // Eliminate
            S pivot = A[i][i];
            for (size_t j = i + 1; j < dim; ++j) {
                S factor = A[j][i] / pivot;
                A[j] -= factor * A[i];
                x[j] -= factor * x[i];
            }
        }

        // Back substitution
        V result(dim);
        for (int i = dim - 1; i >= 0; --i) {
            S sum = x[i];
            for (size_t j = i + 1; j < dim; ++j)
                sum -= A[i][j] * result[j];
            result[i] = sum / A[i][i];
        }

        return result;
    }

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
            for (int i=0; i < ATb.get_dim(); ++i) 
                ATA[i][i] += 1e-10;
            ATb += b[n]*A[n];
        }
        return lin_solve(ATA, ATb);
    }
    
}

#endif
