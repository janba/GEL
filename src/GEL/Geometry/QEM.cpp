/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>

#include "QEM.h"
#include "../CGLA/eigensolution.h"

using namespace CGLA;
using namespace std;

namespace Geometry
{
	Vec3d QEM::opt_pos(double sv_thresh, const CGLA::Vec3d& p0) const
	{
        // Compute eigensolution of the symmetric matrix A. This
        // allows us to factorize it into A = U L U^T and compute
        // the pseudoinverse.
        Mat3x3d U(0),L(0);
        int n = power_eigensolution(A, U, L);
        
        // Unfortunately, eigendecomposition does not find the basis
        // vectors of the 0-space, so we compute either one or two
        // vectors below that span the 0-space.
        switch(n) {
            case 0:
                return p0;
            case 1:
                orthogonal(U[0], U[1], U[2]);
                break;
            case 2:
                U[2] = cross(U[0],U[1]);
                break;
        }
        
        // For each eigenvalue, we compute the corresponding component
        // of either the least squares or least norm solution.
        double limit = abs(sv_thresh * L[0][0]);
        Vec3d x(0);
        for(int i=0;i<3;++i) {
            if(abs(L[i][i])<limit)
                x += U[i] * dot(U[i], p0);
            else
                x -= U[i] * dot(U[i], 0.5*b)/L[i][i];
        }
        return x;
	}
}
