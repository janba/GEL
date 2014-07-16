/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "eigensolution.h"

#include "Mat2x2f.h"
#include "Mat3x3f.h"
#include "Mat4x4f.h"
#include "Mat2x2d.h"
#include "Mat3x3d.h"
#include "Mat4x4d.h"

#include <iostream>

using namespace std; 


namespace
{
    const unsigned int KMAX = 1000000;
    const double EV_THRESH = 1e-6;
}

namespace CGLA
{
    template <class MT>
    int power_eigensolution(const MT& Ap, MT& Q, MT& L, unsigned int max_sol)
    {
        L = MT(0);
        
        typedef typename MT::VectorType VT;
        MT A = Ap;
        unsigned int n = min(MT::get_v_dim(), max_sol);
        
        gel_srand(0);
        for(unsigned int i=0;i<n;++i)
        {
            // Seed the eigenvector estimate
            VT q;
            for (size_t j=0; j<MT::get_v_dim(); ++j) 
                q[j] = gel_rand()/static_cast<double>(GEL_RAND_MAX);
            
            q.normalize();
            double l=123,l_old;
            
            // As long as we haven't reached the max iterations and the
            // eigenvalue has not converged, do
            unsigned int k=0;
            do
            {
                const VT z = A * q;
                double z_len = length(z);
                
                if(z_len < EV_THRESH) return i;
                
                l_old = l;
                l = dot(q, z)>0 ? z_len : -z_len;
                q = z/z_len;
                
                if(++k==KMAX)
                    return i;
            }
            while((fabs(l-l_old) > fabs(EV_THRESH * l)) || k<2);
            
            // Update the solution by adding the eigenvector to Q and
            // the eigenvalue to the diagonal of L.
            Q[i] = q;
            L[i][i] = l;
            
            // Update A by subtracting the subspace represented by the 
            // eigensolution just found. This is called the method of 
            // deflation.
            MT B;
            outer_product(q,q,B);
            A = A - l * B;
        }
        return n;
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
}

