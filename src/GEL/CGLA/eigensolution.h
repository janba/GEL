/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file eigensolution.h
 * @brief Compute eigensolutions for symmetric CGLA matrices.
 */

#ifndef __CGLA_EIGENSOLUTION_H__
#define __CGLA_EIGENSOLUTION_H__

namespace CGLA
{
	/** Use the power method to obtain an eigensolution.
			Given a matrix A, the function returns the number
			of eigensolutions found, and the eigenvectors are
			stored in Q as the rows, and the corresponding 
			values are stored in the diagonal of L upon return of 
			the function.

			The so called power method is used to find the dominant 
			eigenvalue, and the method of deflation is used to find
			the following values. This restricts this function to
			work only on symmetric matrices. 

			DO NOT CALL THIS FUNCTION WITH AN UNSYMMETRIC MATRIX.

			The final argument is the number of solutions to find. If only
			a number of solutions are interesting, use this argument to save
			cycles.
	*/
	template <class MT>
	int power_eigensolution(const MT& A, MT& Q, MT& L, unsigned int max_sol=1000);
}
#endif
