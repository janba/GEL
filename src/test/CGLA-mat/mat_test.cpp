#include <iostream>
#include <algorithm>

#include <GEL/CGLA/Mat4x4d.h>
#include <GEL/CGLA/Mat4x4f.h>
#include <GEL/CGLA/Mat2x2f.h>
#include <GEL/CGLA/Mat3x3f.h>
#include <GEL/CGLA/Mat2x3f.h>
#include <GEL/CGLA/eigensolution.h>

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec2i.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3f.h>

#include <GEL/LinAlg/LapackFunc.h>

using namespace std;
using namespace CGLA;
using namespace LinAlg;

double percentdiff(double a, double b)
{
	return 100.0*fabs(a-b)/fabs(a);
}

void power_eigensolution_test()
{
	gel_srand(0);
	double avg_trace_diff=0;
	double avg_det_diff=0;
	for (int k=0;k<1000;++k)
	{
		Mat4x4d M;
		for(int i=0;i<4;++i) 
			for(int j=i;j<4;++j)
				M[i][j] = M[j][i] = 1e-5*(gel_rand()-GEL_RAND_MAX/2);
		
		Mat4x4d Q,L;
		power_eigensolution(M,Q,L);
		
		double det_diff = percentdiff(determinant(M), L[0][0]*L[1][1]*L[2][2]*L[3][3]) ;
		double trace_diff= percentdiff(trace(M), L[0][0]+L[1][1]+L[2][2]+L[3][3]);
		
		if(det_diff > 1)
		{
			cout << "Matrix " << M << " has Q matrix " <<  Q << " and L " << L << endl;
			CMatrix CM(M);
			CVector CV = SVD(CM);
			cout << "Abs. eigenvalues computed using SVD "; 
			for(int i=0;i<4;++i)
				cout << (CV[i]) <<  " " ;
			cout << endl;
			cout << "Det(M) = " << determinant(M) << " trace(m) = " << trace(M) <<endl;
			
			cout << "Det(M) vs product of eigenvalues, difference = " 
				<< det_diff << " %" << endl;
			avg_det_diff += det_diff;
			cout << "Trace of M vs sum of eigenvalues, difference = " 
				<< trace_diff << " %" << endl;
			avg_trace_diff += trace_diff;
			cout << endl << "-----------------------------------" << endl;
		}
	}
	cout << "average trace diff = " << (avg_trace_diff / 1000) << " %" << endl;
	cout << "average det diff = " << (avg_det_diff / 1000) << " %" << endl;

}


/* This is a non-exhaustive test program for CGLA */
int main()
{

	Mat2x3f m23(Vec3f(1,1,1),Vec3f(2,2,2));
	Mat3x2f m32;
	
	// Try transpose of non-sq matrix
	transpose(m23, m32);
	cout << "2 by 3 " << m23 << endl;
	cout << "Transpose " << m32 << endl;
	
	
	m32[0] = Vec2f(1,2);
	m32[1] = Vec2f(10,20);
	m32[2] = Vec2f(100,200);

	Mat2x2f m22;
	
	// Multiply 2x3 and 3x2 matrices
	mul(m23, m32, m22);
	cout << "Multiplication of non. sq. mat " << m22 << endl;

	// Check matrix mul using * operator for sq matrices.
	Mat2x2f m22_2(1,2,3,4);
	cout << "multiplication of sq. mats " << m22_2 * m22 << endl;
	cout << "Trace " << trace(m22) << endl;
	
	m22_2 += m22 + m22;
	m22_2 -= 2* m22;
	
	Mat4x4f m44(0.0f);
	m44[0][0] = 1;
	m44[1][2] = 2;
	m44[2][1] = 3;
	m44[3][3] = 4;

	Mat4x4f m44_2 = transpose(m44);
	cout << m44 << m44_2 << endl; 

	// Compile should fail if line below is uncommented
	// mul(m23, m44, m22);
	// Compile should fail if line below is uncommented
	// transpose(m23,m44);

	cout << "Determinant of 4x4 " << m44 << determinant(m44) << endl;
	cout << "Determinant of 2x2 " << m22 << determinant(m22) << endl;

	Mat4x4f mna(Vec4f(10,120,10,40),
							Vec4f(43,10,31254,10),
							Vec4f(43,11,54,10),
							Vec4f(0,0,0,1));

	cout << fixed << endl;
	try
		{
			cout << "Invert: " << invert_affine(mna) << endl;
			cout << "Test: " << invert_affine(mna)*mna << endl;
		}
	catch(const Mat4x4fException& me)
		{
			me.print(cout);
		}
	try
		{	
			cout << "Invert " << invert(mna) << endl;
			cout << "test " << invert(mna)*mna << endl;
		}
	catch(const Mat4x4fException& me)
		{
			me.print(cout);
		}

	Vec2f v2(2,3);
	Vec3f v3(1,2,3);
	Mat2x3f m23_2;
	outer_product(v2,v3,m23_2);
	cout << "outer product of " << v2 << " , " << v3 << " = " << m23_2 << endl;
	
	m44 = identity_Mat4x4f();

	Mat3x3f m33(Vec3f(1,2,3), Vec3f(4,5,6), Vec3f(7,8,9));

	copy_matrix(m33, m44);

	cout << "The matrix " << m33 << " is copied to " << m44 << endl;

	cout << "Determinant of singular 4x4 matrix" << endl;

 	Mat4x4d m_odd(Vec4d(0, 373, 139129, 1),
 								Vec4d(-373, -125, 154754, 1),
 								Vec4d(-50.6168, 21.2469, 3013.5, 100),
 								Vec4d(-50.6168, 21.2469, 3013.5, 100));

	cout << m_odd << endl;
	cout << "Determinant: " << determinant(m_odd) << endl;

	Mat4x4d m_odd_t = transpose(m_odd);

	cout << m_odd_t << endl;
	cout << "Det of transpose:" << determinant(m_odd_t) << endl;

	cout << "Determinant of random matrix (MatLab says 0.0491)" << endl;
	

	Mat4x4d Mrand(Vec4d(0.9169,0.3529,0.2028,0.1988),
								Vec4d(0.4103,0.8132,0.1987,0.0153),
								Vec4d(0.8936,0.0099,0.6038,0.7468),
								Vec4d(0.0579,0.1389,0.2722,0.4451));
	
	cout << "det " << determinant(Mrand) << endl;
	cout << "det of transpose " << determinant(transpose(Mrand)) << endl;
	
	cout << "Grand test of power eigensolution scheme " << endl;
	power_eigensolution_test();

	return 0;
}
