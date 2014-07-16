/* LinAlg EigenSolutionSym test program. Andreas Bærentzen 2008. 

	This simple program generates 30 moderately sized dense symmetric real matrices and computes
	the full eigensolution (eigenvalues and eigenvectors). For each matrix and each eigenvector it
	multiplies that eigenvector back onto the matrix and checks the difference between the new vector
	and the eigenvector times the eigenvalue. It also computes the difference between the eigenvalue
	returned and the eigenvalue computed as the norm of a unit length eigenvector multiplied onto the
	matrix.
*/

//disable the stl warning of debug information conflicting with names of longer than 255 char.
#pragma warning( disable :4786 )  

#include <iostream>
#include <vector>
#include <GEL/LinAlg/Matrix.h>
#include <GEL/LinAlg/Vector.h>
#include <GEL/LinAlg/LapackFunc.h>

using namespace std;
using namespace LinAlg;
using namespace CGLA;

int main()
{
		gel_srand(0);
		double max_err_eigen=0;
		double max_err_eigen_vec=0;
		for(int i=0;i<30;++i)
		{
			int N = gel_rand()%50+2;
			CMatrix mat(N,N,0);
			for(int i=0;i<N;++i)
				for(int j=0;j<N;++j)
				{
					double x = gel_rand()/double(GEL_RAND_MAX)-0.5;
					mat.set(i,j,x);
					mat.set(j,i,x);
				}
			CVector b;
			CMatrix Q=mat;
			int info = EigenSolutionsSym(Q,b);
			cout << "Info = " << info << endl;
			for(int i=0;i<b.Length(); ++i)
			{
				CVector q;
				Q.GetRow(q, i);
				q /= q.Norm();

				CVector q2;
				q2 = mat * q;
				double eigen_check = q2.Norm();
				double err_eigen =   (fabs(b[i])-fabs(eigen_check));
				double err_eigen_vec = (q*b[i]-q2).Norm();
				cout << "l(i)-Norm(Q*e(i)) " << err_eigen << endl;
				cout << "Norm(e(i)*l(i)-Q*e(i)) " << err_eigen_vec << endl;
				max_err_eigen = max(err_eigen, max_err_eigen);
				max_err_eigen_vec = max(err_eigen_vec, max_err_eigen_vec);
			}
			cout << "-------" << endl;
		}
		cout << "Maximum eigenvalue error :" << max_err_eigen << endl;
		cout << "Maximum eigenvector error :" << max_err_eigen_vec << endl;


}
