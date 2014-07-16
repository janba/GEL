//disable the stl warning of debug information conflicting with names of longer than 255 char.
#pragma warning( disable :4786 )  

#include <iostream>
#include <vector>
#include <GEL/CGLA/CGLA.h>
#include <GEL/LinAlg/Matrix.h>
#include <GEL/LinAlg/Vector.h>
#include <GEL/LinAlg/LapackFunc.h>

using namespace std;
using namespace LinAlg;
using namespace CGLA;

int main()
{
		// Matrix row and column dimensions
		int N=4,M=4;

		// Fill the matrix with random values.
		gel_srand(0);
		CMatrix mat(N,M,0);
		for(int i=0;i<N;++i)
				for(int j=0;j<M;++j)
						mat.set(i,j,gel_rand()/double(GEL_RAND_MAX)-0.5);

		// Create a vector of random values
		CVector b(N);
		for(int i=0;i<N;++i)
				b.set(i, gel_rand()/double(GEL_RAND_MAX)-0.5);
		
		// Some output
		cout << "Matrix rows, cols = "; 
		cout << mat.Rows() << " , " << mat.Cols() << endl;
		cout << "Matrix: " << mat << endl;
		cout << "Vector: " << b << endl;
		
		// Use lapack function to do least squares solution
		CVector x(M);
		x = LinearLSSolve(mat,b);
		cout << "Least squares solution: " << x << endl;
		
		// Of course if the system has full rank, we can use
		//CVector LinearSolve(const CMatrix& A,const CVector&b);

		// Manual (and less efficient least squares solution ...
		CMatrix mat_T = mat.Transposed(); 
		CMatrix mat2 = mat_T * mat;
		Invert(mat2);
		mat2 = mat2 * mat_T;

		CVector x2(M);
		x2 = mat2 * b;
		cout << "Least squares solution: " << x2 << endl;

}
