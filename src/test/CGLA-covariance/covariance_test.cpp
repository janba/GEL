#include <iostream>

#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Mat3x3f.h>
#include <GEL/CGLA/statistics.h>
#include <GEL/CGLA/eigensolution.h>

using namespace std;
using namespace CGLA;


float frand() { return static_cast<float>(gel_rand())/GEL_RAND_MAX - 0.5f;}

/* This is a non-exhaustive test program for CGLA */
int main()
{
    Vec3f norm(1,0,0);
    norm.normalize();
    
    Vec3f a,b;
    orthogonal(norm,a,b);
    
    cout << "The frame " << endl;
    cout << norm << a << b << endl;
    
    vector<Vec3f> vec;
    gel_srand(0);
    for(int i=0;i<100000;++i)
    {
	Vec3f p = 1.0f*frand()*norm + 2.02f*frand()*a + 3.05f*frand()*b;
	vec.push_back(p);
    }
    
    Mat3x3f A(0);
    Vec3f m = covariance(vec, A);
    
    cout << "Mean and covariance " << endl;
    cout << m << "\n" << A << endl;
    
    
    Mat3x3f Q,L;
    int n = power_eigensolution(A, Q, L);
    
    cout << "The " << n << " eigensolutions are ";
    cout << Q << L << endl;
    
    cout << "Dot products " 
	<< dot(Q[0], Q[1]) << " "
	<< dot(Q[0], Q[2]) << " " 
	<< dot(Q[1], Q[2]) << endl;
    
    cout << "Now with identity" << endl;
    A = identity_Mat3x3f();
    n = power_eigensolution(A, Q, L);
    cout << "The " << n << " eigensolutions are ";
    cout << Q << L << endl;
    
    cout << "Dot products " 
	<< dot(Q[0], Q[1]) << " "
	<< dot(Q[0], Q[2]) << " " 
	<< dot(Q[1], Q[2]) << endl;

	return 0;
}
