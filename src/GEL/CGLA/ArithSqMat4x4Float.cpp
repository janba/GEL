/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ArithSqMat4x4Float.h"
#include "Mat4x4f.h"
#include "Mat4x4d.h"

namespace CGLA {
	
namespace
{
	
	/* Aux func. computes 3x3 determinants. */
	template<class T>
	inline T d3x3f( T a0, T a1, T a2, 
									T b0, T b1, T b2, 
									T c0, T c1, T c2 )
	{
		return a0*b1*c2 +
			a1*b2*c0 +
			a2*b0*c1 -
			a2*b1*c0 -
			a0*b2*c1 -
			a1*b0*c2 ;
	}
}



	template<class V, class M>
	M invert_affine(const ArithSqMat4x4Float<V,M>& this_mat)
	{
		//   The following com[3]e has been copied from a gem in Graphics Gems II by
		//   Kevin Wu.                      
        //   From the EULA: "Using the code is permitted in any program, product, or 
        //   library, non-commercial or commercial. Giving credit is not required, though is a nice gesture"
		//   The function is very fast, but it can only invert affine matrices. An
		//   exception NotAffine is thrown if the matrix is not affine, and another
		//   exception Singular is thrown if the matrix is singular.
		
		typedef typename M::ScalarType ScalarType;

		M new_mat;
		ScalarType det_1;
		ScalarType pos, neg, temp;
      

		if (!(is_tiny(this_mat[3][0]) && 
					is_tiny(this_mat[3][1]) && 
					is_tiny(this_mat[3][2]) && 
					is_tiny(this_mat[3][3]-1.0)))
			throw(Mat4x4fNotAffine("Can only invert affine matrices"));
    
#define ACCUMULATE if (temp >= 0.0) pos += temp; else neg += temp
  
		/*
		 * Calculate the determinant of submatrix A and determine if the
		 * the matrix is singular as limited by the float precision
		 * floating-point this_mat representation.
		 */
  
		pos = neg = 0.0;
		temp =  this_mat[0][0] * this_mat[1][1] * this_mat[2][2];
		ACCUMULATE;
		temp =  this_mat[1][0] * this_mat[2][1] * this_mat[0][2];
		ACCUMULATE;
		temp =  this_mat[2][0] * this_mat[0][1] * this_mat[1][2];
		ACCUMULATE;
		temp = -this_mat[2][0] * this_mat[1][1] * this_mat[0][2];
		ACCUMULATE;
		temp = -this_mat[1][0] * this_mat[0][1] * this_mat[2][2];
		ACCUMULATE;
		temp = -this_mat[0][0] * this_mat[2][1] * this_mat[1][2];
		ACCUMULATE;
		det_1 = pos + neg;
  
		/* Is the submatrix A singular? */
		if ((det_1 == 0.0) || (fabs(det_1 / (pos - neg)) < MINUTE)) 
			{
				/* Mat4x4f M has no inverse */
				throw(Mat4x4fSingular("Tried to invert Singular matrix"));
			}

		else {

			/* Calculate inverse(A) = adj(A) / det(A) */
			det_1 = 1.0 / det_1;
			new_mat[0][0] =   ( this_mat[1][1] * this_mat[2][2] -
													this_mat[2][1] * this_mat[1][2] )
				* det_1;
			new_mat[0][1] = - ( this_mat[0][1] * this_mat[2][2] -
													this_mat[2][1] * this_mat[0][2] )
				* det_1;
			new_mat[0][2] =   ( this_mat[0][1] * this_mat[1][2] -
													this_mat[1][1] * this_mat[0][2] )
				* det_1;
			new_mat[1][0] = - ( this_mat[1][0] * this_mat[2][2] -
													this_mat[2][0] * this_mat[1][2] )
				* det_1;
			new_mat[1][1] =   ( this_mat[0][0] * this_mat[2][2] -
													this_mat[2][0] * this_mat[0][2] )
				* det_1;
			new_mat[1][2] = - ( this_mat[0][0] * this_mat[1][2] -
													this_mat[1][0] * this_mat[0][2] )
				* det_1;
			new_mat[2][0] =   ( this_mat[1][0] * this_mat[2][1] -
													this_mat[2][0] * this_mat[1][1] )
				* det_1;
			new_mat[2][1] = - ( this_mat[0][0] * this_mat[2][1] -
													this_mat[2][0] * this_mat[0][1] )
				* det_1;
			new_mat[2][2] =   ( this_mat[0][0] * this_mat[1][1] -
													this_mat[1][0] * this_mat[0][1] )
				* det_1;

			/* Calculate -C * inverse(A) */
			new_mat[0][3] = - ( this_mat[0][3] * new_mat[0][0] +
													this_mat[1][3] * new_mat[0][1] +
													this_mat[2][3] * new_mat[0][2] );
			new_mat[1][3] = - ( this_mat[0][3] * new_mat[1][0] +
													this_mat[1][3] * new_mat[1][1] +
													this_mat[2][3] * new_mat[1][2] );
			new_mat[2][3] = - ( this_mat[0][3] * new_mat[2][0] +
													this_mat[1][3] * new_mat[2][1] +
													this_mat[2][3] * new_mat[2][2] );
    
			/* Fill in last column */
			new_mat[3][0] = new_mat[3][1] = new_mat[3][2] = 0.0;
			new_mat[3][3] = 1.0;

			return new_mat;
	
		}

#undef ACCUMULATE
	}

	template<class V, class M>
	M adjoint(const ArithSqMat4x4Float<V,M>& in)
	{
		double a1, a2, a3, a4, b1, b2, b3, b4;
		double c1, c2, c3, c4, d1, d2, d3, d4;

		/* assign to individual variable names to aid  */
		/* selecting correct values  */
	
		a1 = in[0][0]; b1 = in[0][1]; 
		c1 = in[0][2]; d1 = in[0][3];

		a2 = in[1][0]; b2 = in[1][1]; 
		c2 = in[1][2]; d2 = in[1][3];

		a3 = in[2][0]; b3 = in[2][1];
		c3 = in[2][2]; d3 = in[2][3];

		a4 = in[3][0]; b4 = in[3][1]; 
		c4 = in[3][2]; d4 = in[3][3];


		/* row column labeling reversed since we transpose rows & columns */
	
		M out;
		out[0][0]  =   d3x3f( b2, b3, b4, c2, c3, c4, d2, d3, d4);
		out[1][0]  = - d3x3f( a2, a3, a4, c2, c3, c4, d2, d3, d4);
		out[2][0]  =   d3x3f( a2, a3, a4, b2, b3, b4, d2, d3, d4);
		out[3][0]  = - d3x3f( a2, a3, a4, b2, b3, b4, c2, c3, c4);
	
		out[0][1]  = - d3x3f( b1, b3, b4, c1, c3, c4, d1, d3, d4);
		out[1][1]  =   d3x3f( a1, a3, a4, c1, c3, c4, d1, d3, d4);
		out[2][1]  = - d3x3f( a1, a3, a4, b1, b3, b4, d1, d3, d4);
		out[3][1]  =   d3x3f( a1, a3, a4, b1, b3, b4, c1, c3, c4);
        
		out[0][2]  =   d3x3f( b1, b2, b4, c1, c2, c4, d1, d2, d4);
		out[1][2]  = - d3x3f( a1, a2, a4, c1, c2, c4, d1, d2, d4);
		out[2][2]  =   d3x3f( a1, a2, a4, b1, b2, b4, d1, d2, d4);
		out[3][2]  = - d3x3f( a1, a2, a4, b1, b2, b4, c1, c2, c4);
	
		out[0][3]  = - d3x3f( b1, b2, b3, c1, c2, c3, d1, d2, d3);
		out[1][3]  =   d3x3f( a1, a2, a3, c1, c2, c3, d1, d2, d3);
		out[2][3]  = - d3x3f( a1, a2, a3, b1, b2, b3, d1, d2, d3);
		out[3][3]  =   d3x3f( a1, a2, a3, b1, b2, b3, c1, c2, c3);

		return out;
	}


	template<class V, class M>
	M invert(const ArithSqMat4x4Float<V,M>& in)
	{
		double det = determinant( in );
		if (is_tiny(det)) 
			throw(Mat4x4fSingular("Tried to invert Singular matrix"));
		M out = adjoint(in);
		out/=det;
		return out;
	}


	template class ArithSqMat4x4Float<Vec4f, Mat4x4f>;
	template Mat4x4f adjoint(const ArithSqMat4x4Float<Vec4f,Mat4x4f>&);
	template Mat4x4f invert(const ArithSqMat4x4Float<Vec4f,Mat4x4f>&);
	template Mat4x4f invert_affine(const ArithSqMat4x4Float<Vec4f,Mat4x4f>&);

	template class ArithSqMat4x4Float<Vec4d, Mat4x4d>;
	template Mat4x4d adjoint(const ArithSqMat4x4Float<Vec4d,Mat4x4d>&);
	template Mat4x4d invert(const ArithSqMat4x4Float<Vec4d,Mat4x4d>&);
	template Mat4x4d invert_affine(const ArithSqMat4x4Float<Vec4d,Mat4x4d>&);


}
