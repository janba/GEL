/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include "LapackFunc.h"

extern "C" {
	extern void dgetrf_(const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info);
	extern void dgetri_(const int *n, double *a, const int *lda, int *ipiv,double *work, int *lwork, int *info);
	extern int dgesvd_(const char *jobu, const char *jobvt, const int *m, const int *n, double *a, const int *lda, double *s,
		double *u, const int *ldu, double *vt, const int *ldvt, double *work, const int *lwork, int *info);
	extern void dgelss_(const int *, const int *, const int *, double *, const int *, double *, const int *, double *, double *, int *, double *, const int *, int *);
	extern void dgesv_(const int *N, const int *NRHS, double *A, const int *LDA, int *IPIV, double *B, const int *LDB, int *INFO);
	extern void dposv_(const char *UPLO, const int *N, const int *NRHS, double *A, const int *LDA, double *B, const int *LDB, int *INFO);
	extern void dgeqrf_(const int *m, const int *n, double* A,const int* lda,double*tau,double *work, const int* lwork, int* info);
	extern void dorgqr_(const int *m, const int *n, const int *k, double* A, const int *lda,const double *tau,double* work, const int *lwork, int *info);
	extern void dgerqf_(const int *m, const int *n, double* A,const int* lda,double*tau,double *work, const int* lwork, int* info);
	extern void dorgrq_(const int *m, const int *n, const int *k, double* A, const int *lda,const double *tau,double* work, const int *lwork, int *info);

	int dsysv_(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, double *work, int *lwork, int *info);
	
	int dsytrd_(char *uplo, int *n, double *a, int * lda, double *d__, double *e, double *tau, double *work, int *lwork, int *info);
	int dsteqr_(char *compz, int *n, double *d__, double *e, double *z__, int *ldz, double *work, int *info);
    int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

}

namespace{
	
	template <class T>
		inline T MIN(const T& a,const T& b)
	{
		return a<b?a:b;
	}
	
	template <class T>
		inline T MAX(const T& a,const T& b)
	{
		return a>b?a:b;
	}
	
}

namespace LinAlg
{

void SVD(const CMatrix& A,
		 CMatrix& U,
		 CVector& s,
		 CMatrix& V)
{
	//Notice, that the Lapack/Fortran CMatrix representation is the transposed version of the C/C++
	CMatrix a;
	A.Transposed(a);
	
	int nRows=A.Rows();
	int nCols=A.Cols();
	int info;
	int lda = MAX(nRows,1);
	int lwork= 6*MIN(nRows,nCols) + nRows + nCols;
	
	s.Resize(MIN(nRows,nCols));
	double *work = new double[lwork];
	
	char jobu='A';
	int ldu=nRows;
	U.Resize(nRows,nRows);
	
	char jobvt='A';
	int ldvt=nCols;
	V.Resize(nCols,nCols);
	
	
	dgesvd_ (&jobu, &jobvt, &nRows, &nCols, a[0], &lda, &(s[0]), U[0], &ldu, V[0], &ldvt,work, &lwork, &info);
	
	assert(info>=0);
	
	delete [] work;
	U.Transpose();
}

void SVD(const CMatrix& A,
		 CMatrix& U,
		 CMatrix& S,
		 CMatrix& V)
{
	CVector s;
	SVD(A,U,s,V);
	S.Resize(A.Rows(),A.Cols());
	S=0;
	for(int cS=MIN(A.Rows(),A.Cols())-1;cS>=0;cS--)
	{
		S[cS][cS]=s[cS];
	}
}

CVector SVD(const CMatrix& A)
{
	CMatrix a;
	A.Transposed(a);
	
	int nRows=A.Rows();
	int nCols=A.Cols();
	int info;
	int lda = MAX(nRows,1);
	int lwork= 6*MIN(nRows,nCols) + nRows + nCols;
	
	CVector s(MIN(nRows,nCols));
	double *work = new double[lwork];
	
	char jobu='N';
	int ldu=1;
	
	char jobvt='N';
	int ldvt=1;
	
	dgesvd_ (&jobu, &jobvt, &nRows, &nCols, a[0], &lda, &(s[0]), NULL, &ldu, NULL, &ldvt,work, &lwork, &info);
	
	assert(info>=0);
	
	delete [] work;
	
	return s;
}


void LinearSolve(const CMatrix& A,
				 const CVector&b,
				 CVector& x)
{
	assert(A.Rows()==b.Length());
	assert(A.Rows()==A.Cols());
	
	CMatrix a;
	A.Transposed(a);
	x=b;
	int nRows=A.Rows();
	int nrhs=1;	//only one CVector, change here to make b a CMatrix.
	int info;
	int *ipiv =new int[nRows];
	
	dgesv_(&nRows, &nrhs, a[0], &nRows, ipiv, &(x[0]), &nRows, &info);
	
	assert(info==0);
	delete [] ipiv;
}

CVector LinearSolve(const CMatrix& A,
				   const CVector&b)
{
	CVector x;
	LinearSolve(A,b,x);
	return x;
}


void LinearSolveSPD(const CMatrix& A,
					const CVector&b,
					CVector& x)
{
	assert(A.Rows()==b.Length());
	assert(A.Rows()==A.Cols());
	
	CMatrix a(A);
	x=b;
	
	char uplo='U';
	int nRows=A.Rows();
	int info;
	int nrhs=1;
	
	dposv_(&uplo, &nRows, &nrhs, a[0], &nRows, &(x[0]), &nRows, &info);
	
	assert(info==0);
}

void LinearSolveSym(const CMatrix& A,
										const CVector&b,
										CVector& x)
{
	assert(A.Rows()==b.Length());
	assert(A.Rows()==A.Cols());
	
	CMatrix a(A);
	x=b;
	
	char uplo='U';
	int nRows=A.Rows();
	int nCols=A.Cols();
	int info;
	int nrhs=1;
	
	int lwork= 10*6*MIN(nRows,nCols) + nRows + nCols;
	double *work = new double[lwork];

	std::vector<int> ipiv(A.Rows());
	dsysv_(&uplo, &nRows, &nrhs, a[0], &nRows, &ipiv[0], &(x[0]), &nRows, 
				 work, &lwork, &info);
	delete [] work;
	assert(info==0);
}


CVector LinearSolveSPD(const CMatrix& A,const CVector&b)
{
	CVector x;
	LinearSolveSPD(A,b,x);
	return x;
}

void LinearLSSolve(const CMatrix& A,
				   const CVector& b,
				   CVector& x)
{
	assert(A.Rows()==b.Length());
	
	int nRows=A.Rows();
	int nCols=A.Cols();
	CMatrix a;
	A.Transposed(a);
	
	int ldb=MAX(nRows,nCols);
	double* BX=new double[ldb];
	memcpy(BX,&(b[0]),nRows*sizeof(double));
	
	int nrhs=1;
	double* s=new double[MIN(nRows,nCols)];
	double rcond = -1.0; // using machine precision
	int info,rank;
	int lwork = 5*nRows*nCols + 1; // larger than necessary
	double *work = new double[lwork];
	
	dgelss_(&nRows, &nCols, &nrhs, a[0], &nRows, BX, &ldb, s, &rcond, &rank, work, &lwork, &info);
	
	x.Resize(nCols);
	memcpy(&(x[0]),BX,nCols*sizeof(double));
	
	
	assert(info==0);
	
    delete [] BX;
	delete [] s;
	delete [] work;
}

CVector LinearLSSolve(const CMatrix& A,
					 const CVector& b)
{
	CVector x;
	LinearLSSolve(A,b,x);
	return x;
}

void Invert(CMatrix& A)
{
	assert(A.Rows()==A.Cols());
	int nRows=A.Rows();
	int info;
	int *ipiv = new int[nRows];
	//Perform th LU factorization of A
	dgetrf_ (&nRows, &nRows, A[0], &nRows, ipiv, &info);
	if (info != 0) { //An error occured
		delete [] ipiv;
		assert(info==0);	//info will be < 0 if A is rank deficient.
		return;
	}

	//Calculate the inverse
	int lwork = nRows * 5; 
	double *work = new double[lwork];
	dgetri_ (&nRows, A[0], &nRows, ipiv, work, &lwork, &info);
	delete [] ipiv;
	delete [] work;
	assert(info==0);
}

void Inverted(const CMatrix& A,
			  CMatrix& InvA)
{
	InvA=A;
	Invert(InvA);
}

CMatrix Inverted(const CMatrix& A)
{
	CMatrix InvA(A);
	Invert(InvA);
	return InvA;
}

void QRfact(const CMatrix& A,CMatrix& Q, CMatrix& R)
{

	
	int nRows=A.Rows();
	int nCols=A.Cols();
	assert(A.Rows()>0 && A.Cols()>0); //herby lda=nRows.
	int lwork= 6*nCols;
	int info;
	int k=MIN(nRows,nCols);


	A.Transposed(R);

	double *work = new double[lwork];


	CVector TAU(k);

	dgeqrf_(&nRows,&nCols,R[0],&nRows,&TAU[0],work,&lwork,&info);
	assert(info==0);


	Q=R;

	nCols=MIN(nCols,nRows);

	dorgqr_(&nRows,&nCols,&k,Q[0],&nRows,&TAU[0],work,&lwork,&info);
	assert(info==0);


	delete [] work;

	R.Transpose();
	Q.Transpose();

	//Set the zeros of R
	double * pRow;
	int ColStop;
	for(int cRow=1;cRow<R.Rows();cRow++)
	{
		pRow=R[cRow];
		ColStop=MIN(cRow,R.Cols());
		for(int cCol=0;cCol<ColStop;cCol++)
		{
			pRow[cCol]=0;
		}
	}

}


void RQfact(const CMatrix& A,CMatrix& R, CMatrix& Q)
{
	int nRows=A.Rows();
	int nCols=A.Cols();
	assert(A.Rows()>0 && A.Cols()>0); //herby lda=nRows.
	int lwork= 6*nCols;
	int info;
	int k=MIN(nRows,nCols);


	A.Transposed(R);

	double *work = new double[lwork];


	CVector TAU(k);

	dgerqf_(&nRows,&nCols,R[0],&nRows,&TAU[0],work,&lwork,&info);
	assert(info==0);

	Q=R;

	nCols=MAX(nCols,nRows);

	dorgrq_(&nRows,&nCols,&k,Q[0],&nRows,&TAU[0],work,&lwork,&info);
	assert(info==0);

	
	delete [] work;

	R.Transpose();
	Q.Transpose();

	//Set the zeros of R
	double * pRow;
	int ColStop;
	for(int cRow=1;cRow<R.Rows();cRow++)
	{
		pRow=R[cRow];
		ColStop=MIN(cRow,R.Cols());
		for(int cCol=0;cCol<ColStop;cCol++)
		{
			pRow[cCol]=0;
		}
	}

}

int EigenSolutionsSym(CMatrix& Q, CVector& b)
{
	assert(Q.Rows() == Q.Cols());
	int n = Q.Rows();
	char jobz='V', uplo='U';
	int info, lwork=6*n;
	b.Resize(n);
	double* work = new double[lwork];
	
	dsyev_(&jobz, &uplo, &n, Q[0], &n, &b[0], work, &lwork, &info);
    
    delete [] work;
    
	return info;
}

}
