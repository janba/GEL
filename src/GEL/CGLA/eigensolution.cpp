/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/CGLA/eigensolution.h>

#include <GEL/CGLA/Mat.h>

#include <cmath>
#include <iostream>

using namespace std; 


namespace
{
    const unsigned int KMAX = 1000000;
    const double EV_THRESH = 1e-6;

    // Steps after which a stuck residual means {σ, -σ} rather than a slow gap.
    const unsigned int STALL_STEPS = 32;
}

namespace CGLA
{
    template <class MT>
    int power_eigensolution(const MT& Ap, MT& Q, MT& L, unsigned int max_sol)
    {
        L = MT(0);
        
        typedef typename MT::VectorType VT;
        MT A = Ap;
        const unsigned int dim = MT::get_v_dim();
        unsigned int n = min(dim, max_sol);
        
        gel_srand(0);
        for(unsigned int i=0;i<n;++i)
        {
            // shift == 0 is ordinary power iteration, which finds the largest
            // magnitude first. A pair {σ, -σ} makes ||Aq|| independent of the
            // direction, so the vector oscillates and the residual never
            // falls. The second attempt adds shift*I, whose eigenvalues have
            // different magnitudes, then subtracts shift again. Deflation is
            // always applied to the unshifted matrix.
            double shift = 0;
            VT q;
            double rayleigh = 0;
            bool found = false;

            for (int attempt = 0; attempt < 2; ++attempt) {
                MT M = A;
                if (shift != 0) {
                    for (unsigned int r = 0; r < dim; ++r)
                        M[r][r] += shift;
                }

                for (unsigned int j = 0; j < dim; ++j)
                    q[j] = gel_rand01();
                q.normalize();

                double res0 = 0;
                bool stalled = false;
                found = false;
                for (unsigned int k = 0; k < KMAX; ++k) {
                    const VT z = M * q;
                    const double z_len = length(z);

                    // Nothing of significant magnitude left in the unshifted
                    // matrix. A shifted matrix hides that, so only the first
                    // attempt may decide the spectrum is exhausted.
                    if (z_len < EV_THRESH) {
                        if (shift == 0)
                            return i;
                        break;
                    }

                    rayleigh = dot(q, z);
                    const double res = length(z - rayleigh * q);
                    if (res <= EV_THRESH * (1.0 + fabs(rayleigh))) {
                        found = true;
                        break;
                    }
                    if (k == 0)
                        res0 = res;

                    q = z / z_len;

                    if (attempt == 0 && k == STALL_STEPS && res > 0.99 * res0) {
                        shift = 0.25 * z_len;
                        stalled = true;
                        break;
                    }
                }
                if (found)
                    break;
                if (!stalled)
                    return i;
            }
            if (!found)
                return i;

            const double lambda = rayleigh - shift;
            if (fabs(lambda) < EV_THRESH)
                return i;

            // Update the solution by adding the eigenvector to Q and
            // the eigenvalue to the diagonal of L.
            Q[i] = q;
            L[i][i] = lambda;
            
            // Update A by subtracting the subspace represented by the 
            // eigensolution just found. This is called the method of 
            // deflation.
            MT B;
            outer_product(q,q,B);
            A = A - lambda * B;
        }
        return n;
    }

    template <class MT>
    int eigensolution_2x2(const MT& A, MT& Q, MT& L)
    {
        using V = typename MT::VectorType;
        using S = typename MT::ScalarType;
        static_assert(MT::get_v_dim() == 2, "eigensolution_2x2 requires a 2x2 matrix");

        const S a = A[0][0];
        const S off = S(0.5) * (A[0][1] + A[1][0]);
        const S c = A[1][1];
        const S trace = a + c;
        const S diff = a - c;
        const S disc = std::sqrt(diff * diff + S(4) * off * off);
        const S l_max = S(0.5) * (trace + disc);
        const S l_min = S(0.5) * (trace - disc);

        // (A - λI)q = 0. The row with the larger diagonal residual gives a
        // stable null vector; a pure multiple of the identity has neither.
        const S row0 = std::fabs(a - l_max);
        const S row1 = std::fabs(c - l_max);
        V q0 = row0 >= row1 ? V(off, l_max - a) : V(l_max - c, off);
        if (!(length(q0) > S(0)))
            q0 = V(S(1), S(0));
        else
            q0.normalize();
        const V q1(-q0[1], q0[0]);

        Q = MT(0);
        L = MT(0);
        Q[0] = q0;
        Q[1] = q1;
        L[0][0] = l_max;
        L[1][1] = l_min;
        return 2;
    }
    
    /* There is no reason to put this template in a header file, since 
     we will only use it on matrices defined in CGLA. Instead, we 
     explicitly instantiate the function for the square matrices
     of CGLA */
    template int power_eigensolution<Mat2x2f>(const Mat2x2f&,
                                              Mat2x2f&,Mat2x2f&,unsigned int);
	
    template int power_eigensolution<Mat3x3f>(const Mat3x3f&,
                                              Mat3x3f&,Mat3x3f&,unsigned int);
    template int power_eigensolution<Mat4x4f>(const Mat4x4f&,
                                              Mat4x4f&,Mat4x4f&,unsigned int);
    template int power_eigensolution<Mat2x2d>(const Mat2x2d&,
                                              Mat2x2d&,Mat2x2d&,unsigned int);
	
    template int power_eigensolution<Mat3x3d>(const Mat3x3d&,
                                              Mat3x3d&,Mat3x3d&,unsigned int);
    template int power_eigensolution<Mat4x4d>(const Mat4x4d&,
                                              Mat4x4d&,Mat4x4d&,unsigned int);

    template int eigensolution_2x2<Mat2x2f>(const Mat2x2f&, Mat2x2f&, Mat2x2f&);
    template int eigensolution_2x2<Mat2x2d>(const Mat2x2d&, Mat2x2d&, Mat2x2d&);
}

