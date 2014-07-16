/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ArithSqMat3x3Float.h"
#include "Mat3x3f.h"
#include "Mat3x3d.h"

namespace CGLA {

	template<class V, class M>
	M invert(const ArithSqMat3x3Float<V,M>& _a)
	{
        // From Graphcis Gems IV, Jean-Francois Doue, C++ Vector and Matrix Algebra Routines.
        // From the EULA "Using the code is permitted in any program, product, or library, non-commercial 
        // or commercial. Giving credit is not required, though is a nice gesture"
		M a(_a[0], _a[1], _a[2]);
		M	b;
		b.identity();

		unsigned int  i, j, i1;

		for (j=0; j<3; j++) 
			{   
				i1 = j;                 // Row with largest pivot candidate
				for (i=j+1; i<3; i++)
					if (fabs(a[i][j]) > fabs(a[i1][j]))
						i1 = i;
			
				// Swap rows i1 and j in a and b to put pivot on diagonal
				V a_tmp = a[i1];
				a[i1] = a[j];
				a[j]  = a_tmp;
			
				V b_tmp = b[i1];
				b[i1] = b[j];
				b[j]  = b_tmp;
		
				// Scale row j to have a unit diagonal
				if (a[j][j] == 0.0)
					throw(Mat3x3fSingular("Tried to invert Singular matrix"));
			
				b[j] /= a[j][j];
				a[j] /= a[j][j];
			
				// Eliminate off-diagonal elems in col j of a, doing identical ops to b
				for (i=0; i<3; i++)
					if (i!=j) 
						{
							b[i] -= a[i][j] * b[j];
							a[i] -= a[i][j] * a[j];
						}
			}
		return b;
	}                                                                               
	template class ArithSqMat3x3Float<Vec3f,Mat3x3f>;
	template Mat3x3f invert(const ArithSqMat3x3Float<Vec3f,Mat3x3f>&);

	template class ArithSqMat3x3Float<Vec3d,Mat3x3d>;
	template Mat3x3d invert(const ArithSqMat3x3Float<Vec3d,Mat3x3d>&);
}
