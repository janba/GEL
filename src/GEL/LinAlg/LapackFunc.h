/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef LAPACKFUNC_H_HAA_AGUST_2001
#define LAPACKFUNC_H_HAA_AGUST_2001

//#if defined(_MSC_VER)
//    #if defined(_DEBUG)
//        #pragma message("Note: including lib: lapackd.lib and ignoring defaultlib : LIBC\n")
//        #pragma comment(lib, "lapackd.lib")
//    #else
//        #pragma message("Note: including lib: lapack.lib and ignoring defaultlib : LIBC\n")
//        #pragma comment(lib, "lapack.lib") 
//    #endif
//    #pragma comment(linker, "/NODEFAULTLIB:LIBC.LIB")
//#endif

#include "Matrix.h"
#include "Vector.h"

namespace LinAlg
{

/*!
\file LapackFunc.h
\brief Interface to some of the LAPACK functionality.

These are functions which more or less directly interface with the
Lapack provided algorithms.

For indepth reference to the LAPACK functions see:
  LAPACK Users' Guide - 3rd Edition,
  by E. Anderson et al.,
  ISBN 0-89871-447-8,
  Published by SIAM,

This book is also available at: \URL{http://www.netlib.org/lapack/lug/lapack_lug.html}

The official LAPACK sites where from the source can be downloaded are:
  \URL{http://www.netlib.org/clapack/} and
  \URL{http://www.netlib.org/lapack/}


NB: When running this in MS Visual C++ it is usually required to set the 
multithread "\MD" compiler option. This is to ensure correct linkage to the 
precompiled library "clapack.lib" and/or "clapackDB.lib".


\author  Henrik Aanæs
\version Aug 2001
*/

/*!
\name Singular Value Decomposition SVD 

These functions perform the Singular Value Decomposition SVD of 
the MxN matrix A. The SVD is defined by:

  A=U*S*V^T

where:
- U is a M by M orthogonal matrix 
- V is a N by N orthogonal matrix 
- S is a M by N diaggonal matrix. The values in the diagonal are the singular values


\param  A the matrix to perform SVD on  
\return U will be resized if it is does not have the correct dimensions
\return V will be resized if it is does not have the correct dimensions
\return S will be resized if it is does not have the correct dimensions. 
\exception assert(info==0) for Lapack. Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/  
//@{ 
///SVD of A, where the singular values are returned in a Vector.
void SVD(const CMatrix& A,CMatrix& U,CVector& s,CMatrix& V);
///SVD of A, where the singular values are returned in a 'diagonal' Matrix.
void SVD(const CMatrix& A,CMatrix& U,CMatrix& S,CMatrix& V);
///SVD of A, returning only the singular values in a Vector.
CVector SVD(const CMatrix& A);
//@}


/*!
\name Linear Equations
These functions solve the system of linear equations

  A*x=b 

for x, where:
- A is a N by N matrix 
- b is a N vector
- x is a N vector

There a speceilaized functions for symetric positive definite (SPD) 
matrices yeilding better performance. These are denote by SPD in 
there function name.

\param A the NxN square matrix
\param b the N vector
\return x will be resized if it is does not have the correct dimensions
\exception assert(info==0) for Lapack. Add a throw statement later.
\exception assert(A.Row()==A.Col()). Add a throw statement later.
\exception assert(A.Row()==b.Length()). Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/  
//@{
///Solves Ax=b for x.
void LinearSolve(const CMatrix& A,const CVector&b,CVector& x);
///Solves Ax=b for x and returns x.
CVector LinearSolve(const CMatrix& A,const CVector&b);
///Solves Ax=b for x, where A is SPD.
void LinearSolveSPD(const CMatrix& A,const CVector&b,CVector& x);
///Solves Ax=b for x and returns x, where A is SPD.
CVector LinearSolveSPD(const CMatrix& A,const CVector&b);
//@}

void LinearSolveSym(const CMatrix& A,
										const CVector&b,
										CVector& x);

/**
\name Linear Least Squares
These functions solve the Linear Least Squares problem:

  min_x ||Ax-b||^2

for x, where:
- || || denotes the 2-norm
- A	is a M by N matrix. For a well formed M>=N and rank (A)=N. See below.
- b	is a M vector.
- x	is a N vector 

If the solution is not \em well \em formed the algorithm provided will find a
solution, x, which is not unique, but which sets the objective function
to 0. The reson being that the underlining algorithm works by SVD.

\param A the MxN matrix
\param b the M vector
\return x will be resized if it is does not have the correct dimensions
\exception assert(info==0) for Lapack. Add a throw statement later.
\exception assert(A.Rows()==b.Length());. Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/  
//@{
///Solves the Linear Least Squares problem min_x ||Ax=b||^2 for x.
void LinearLSSolve(const CMatrix& A,const CVector&b,CVector& x);
///Solves the Linear Least Squares problem min_x ||Ax=b||^2 for x, and returnes x.
CVector LinearLSSolve(const CMatrix& A,const CVector&b);
//@}

/**
\name Matrix Inversion
These functions inverts the square matrix A. This matrix A must have
full rank. 
\param A square matrix
\return InvA the invers of A for one instance.
\exception assert(info==0) for Lapack. This wil among others happen if A is rank deficient. Add a throw statement later. 
\exception assert(A.Rows()==A.Cols()). Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/  
//@{
///Invertes the square matrix A. That is here A is altered as opposed to the other Invert functions.
void Invert(CMatrix& A);
/// Returns the inverse of the square matrix A in InvA.
void Inverted(const CMatrix& A,CMatrix& InvA);
/// Returns the inverse of the square matrix A.
CMatrix Inverted(const CMatrix& A);
//@}


/**
\name QR Factorization
This function returns the QR factorization of A, such that Q*R=A where 
Q is a orthonormal matrix and R is an upper triangular matrix. However, 
in the case of A.Col()>A.Row(), the last A.Col-A.Row columns of Q are 
'carbage' and as such not part of a orthonormal matrix.
 \param A  the input matrix
\return Q an orthonormal matrix. (See above)
\return R an upper triangular matrix.
\exception assert(info==0) for Lapack. This wil among others happen if A is rank deficient. Add a throw statement later. 
\exception assert(A.Rows()>0 && A.Cols()>0). Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/ 
//@{ 
void QRfact(const CMatrix& A,CMatrix& Q, CMatrix& R);
//@}


/**
\name RQ Factorization
This function returns the RQ factorization of A, such that R*Q=A where 
Q is a orthonormal matrix and R is an upper triangular matrix. However, 
in the case of A not beeing a square matrix, there might be some fuck up of Q.
 \param A  the input matrix
\return Q an orthonormal matrix. (See above)
\return R an upper triangular matrix.
\exception assert(info==0) for Lapack. This wil among others happen if A is rank deficient. Add a throw statement later. 
\exception assert(A.Rows()>0 && A.Cols()>0). Add a throw statement later.
\version  Aug 2001 
\author  Henrik Aanæs
*/ 
//@{ 
void RQfact(const CMatrix& A,CMatrix& R, CMatrix& Q);
//@}

/**
\name Find eigensolutions of a symmetric real matrix.
This function accepts a real symmetric matrix Q and a vector b.
When the function returns, the eigenvalues of the matrix Q will be stored in b and
the eigenvectors form the rows of Q. This function is based on the
Lapack function dsyev, and returns its info code. A code of 0 indicates
success, and a code < 0 indicates an error. Probably Q is not a real symmetric matrix.
If the code is > 0 "the algorithm failed  to  converge;  code  off-diagonal  elements  
of an intermediate tridiagonal form did not converge to zero." Presumably this means that 
code contains the number of eigenvalues which are ok. 
\author Andreas B¾rentzen.
*/
//@{
int EigenSolutionsSym(CMatrix& Q, CVector& b);
//@}
}

#endif // !defined(LAPACKFUNC_H_HAA_AGUST_2001)
