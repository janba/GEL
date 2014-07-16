#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Quaternion.h>

using namespace std;
using namespace CGLA;

/* CGLA example no. 1:
	 Computing a normal from the three vertices
	 of a triangle */

void example1()
{
	// Define triangle vertices
	Vec3f p0(10,10,10);
	Vec3f p1(20,10,10);
	Vec3f p2(10,20,10);

	// Compute normal
	Vec3f n = normalize(cross(p1-p0,p2-p0));

	// Vectors can be printed using normal facilities
	cout << n << endl;
}

/* CGLA example no. 2:
	 Using spherical coordinates.
*/
void example2()
{
	// Create a null vector
	Vec3f p;
 	p.set_spherical(0.955317f, 3.1415926f/4.0f , 1.0f);
	cout << " p = " << p << ", p.length() = " << p.length() << endl;
}

void example3()
{
	// Construct a quaternion that rotates the x-axis into the 
	// y-axis.
	Quaternion q;
	q.make_rot(Vec3f(1,0,0), Vec3f(0,1,0));
	
	// Create a 4 by 4 Matrix
	Mat4x4f m = translation_Mat4x4f(Vec3f(1,2,3));

	// Multiply the quaternion rotation onto m.
	m *= q.get_Mat4x4f();

	// Create a vector (1,1,1)
	Vec3f p(1);

	// Multiply p onto m as a point. That is we assign 1 to p's w coordinate
	// and then throw away the w coordinate after the transformation.
	Vec3f p2 = m.mul_3D_point(p);

	cout << p2 << endl;
}


void example4()
{
	Vec3f p1(0.8660f,.4f,-.3f);
	Vec3f p2(-0.8660f,.4f,-.3f);
	Vec3f p3(0.0f,.7f,-.7141f);

	Vec3f n1 = cross(p1,p2);
	Vec3f n2 = cross(p2,p3);
	Vec3f n3 = cross(p3,p1);
	Vec3f n  = normalize(n1+n2+n3);
	cout << n1 << n1.length() << endl
			 << n2 << n2.length() << endl
			 << n3 << n3.length() << endl
			 << dot(n,Vec3f(0,0,1)) << endl;
}

int main()
{
// 	example1();
// 	example2();
// 	example3();
	example4();
}
