#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <climits>

#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec2i.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3f.h> 
#include <GEL/CGLA/Vec3d.h> 
#include <GEL/CGLA/Vec3uc.h>
#include <GEL/CGLA/Vec3usi.h>
#include <GEL/CGLA/Vec4f.h>
#include <GEL/CGLA/Vec4d.h>
#include <GEL/CGLA/Vec4uc.h>
#include <GEL/CGLA/Quatf.h>
#include <GEL/CGLA/Quatd.h>
#include <GEL/Util/Timer.h>

using namespace std;
using namespace CGLA;
using namespace Util;

namespace
{
  // To find a pseudo-random number in the interval [0,1[
  double my_random()
  {
    return gel_rand()/(static_cast<double>(GEL_RAND_MAX) + 1.0);
  }
}

/* This is a non-exhaustive test program for CGLA */
int main()
{
    cout << "Note that some of these tests appear to fail although everything works," << endl;
    cout << "e.g. because floating point to double conversion and back does not yield original number" << endl;
  int success = 0;
  Timer t;

  gel_srand(time(0));

  t.start();

  ////////////////////////////////////////////////
  //                V e c 2 d
  ////////////////////////////////////////////////

  // Constructors

  {
  // TODO: The following should be done for all constructors in the remainder of the file
#ifndef NDEBUG
    Vec2d x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]);
#else
    success += 1;
#endif

    double xi1 = my_random(), xi2 = my_random();
    Vec2d x2(xi1, xi2);
    success += x2[0] == xi1 && x2[1] == xi2;

    int xii1 = gel_rand(), xii2 = gel_rand();
    Vec2i i1(xii1, xii2);
    Vec2d x3(i1);
    success += x3[0] == xii1 && x3[1] == xii2;

    float xiii1 = my_random(), xiii2 = my_random();
    Vec2f f1(xiii1, xiii2);
    Vec2d x4(f1);
    success += x4[0] == xiii1 && x4[1] == xiii2;

    double xiiii = my_random();
    Vec2d x5(xiiii);
    success += x5[0] == xiiii && x5[1] == xiiii;
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec2d Constructors";
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec2d x1;    

    success += x1.get_dim() == 2;

    double xi1 = my_random(), xi2 = my_random();
    x1.set(xi1, xi2);
    success += x1[0] == xi1 && x1[1] == xi2;

    success += x1.get() == &x1[0];

    double temp = BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	      temp = x1[i];
    success += temp == x1.min_coord();

    temp = -BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	      temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec2d Data manipulation";
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    double xi1 = my_random(), xi2 = my_random();
    while(xi1 == xi2)
    {
      xi1 = my_random();
      xi2 = my_random();
    }

    Vec2d x1(xi1, xi2), x2(xi1, xi2), x3(xi1), x4(xi2);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);

    if(xi1 < xi2)
    {
      success += x3.all_l(x4);
      success += !(x3.all_l(x2));
      success += x3.all_le(x2);
      success += !(x4.all_le(x3));
      success += x4.all_g(x3);
      success += !(x4.all_g(x2));
      success += x4.all_ge(x2);
      success += !(x3.all_ge(x4));
    }
    else 
    {
      success += x4.all_l(x3);
      success += !(x4.all_l(x2));
      success += x4.all_le(x2);
      success += !(x3.all_le(x4));
      success += x3.all_g(x4);
      success += !(x3.all_g(x2));
      success += x3.all_ge(x2);
      success += !(x4.all_ge(x3));
    }
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec2d Comparison operators";
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    double xi1 = my_random(), xi2 = my_random();
    Vec2d x1(xi1, xi2);
    double xi3 = my_random();
    x1 *= xi3;
    success += x1[0] == xi1*xi3 && x1[1] == xi2*xi3;

    while(xi3 == 0)
      xi3 = my_random();
    x1 = Vec2d(xi1, xi2);
    x1 /= xi3;
    success += abs(x1[0] - xi1/xi3) < 1.0e-15 && abs(x1[1] - xi2/xi3) < 1.0e-15;

    x1 = Vec2d(xi1, xi2);
    x1 += xi3;
    success += x1[0] == xi1 + xi3 && x1[1] == xi2 + xi3;

    x1 = Vec2d(xi1, xi2);
    x1 -= xi3;
    success += x1[0] == xi1 - xi3 && x1[1] == xi2 - xi3;

    double xii1 = my_random(), xii2 = my_random();
    Vec2d x2(xii1, xii2);
    x1 = Vec2d(xi1, xi2);
    x2 *= x1;
    success += x2[0] == xi1*xii1 && x2[1] == xi2*xii2;

    while(xi1 == 0)
      xi1 = my_random();
    while(xi2 == 0)
      xi2 = my_random();
    x1 = Vec2d(xi1, xi2);
    x2 = Vec2d(xii1, xii2);
    x2 /= x1;
    success += abs(x2[0] - xii1/xi1) < 1.0e-15 && abs(x2[1] - xii2/xi2) < 1.0e-15;

    x2 = Vec2d(xii1, xii2);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2;
    
    x2 = Vec2d(xii1, xii2);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec2d Assignment operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Unary operators

  {
    double xi1 = my_random(), xi2 = my_random();
    Vec2d x1 = -Vec2d(xi1, xi2);    
    success += x1[0] == -xi1 && x1[1] == -xi2;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec2d Unary operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Binary operators

  {
    double xi1 = my_random(), xi2 = my_random();
    Vec2d x1(xi1, xi2);
    double xii1 = my_random(), xii2 = my_random();
    while(xii1 == 0)
      xii1 = my_random();
    while(xii2 == 0)
      xii2 = my_random();
    Vec2d x2(xii1, xii2);
    Vec2d x3 = x1*x2;
    success += x3[0] == xi1*xii1 && x3[1] == xi2*xii2;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2;
    
    x3 = x1/x2;
    success += abs(x3[0] - xi1/xii1) < 1.0e-15 && abs(x3[1] - xi2/xii2) < 1.0e-15;

    double xi3 = my_random();
    x3 = x1*xi3;
    success += x3[0] == xi1*xi3 && x3[1] == xi2*xi3;
    
    x3 = xi3*x1;
    success += x3[0] == xi1*xi3 && x3[1] == xi2*xi3;

    float xi4 = my_random();
    x3 = xi4*x1;
    success += x3[0] == xi1*xi4 && x3[1] == xi2*xi4;
    
    x3 = x1*xi4;
    success += x3[0] == xi1*xi4 && x3[1] == xi2*xi4;

    int xi5 = gel_rand();
    x3 = xi5*x1;
    success += x3[0] == xi1*xi5 && x3[1] == xi2*xi5;
    
    x3 = x1*xi5;
    success += x3[0] == xi1*xi5 && x3[1] == xi2*xi5;

    while(xi3 == 0)
      xi3 = my_random();
    x3 = x1/xi3;
    success += abs(x3[0] - xi1/xi3) < 1.0e-15 && abs(x3[1] - xi2/xi3) < 1.0e-15;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec2d Binary operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Vector operations

  {
    double xi1 = my_random(), xi2 = my_random();
    Vec2d x1(xi1, xi2);
    double xii1 = my_random(), xii2 = my_random();
    Vec2d x2(xii1, xii2);
    double x = dot(x1, x2);
    success += x == xi1*xii1 + xi2*xii2;

    x = sqr_length(x1);
    success += x == xi1*xi1 + xi2*xi2;

    Vec2d x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) && x3[1] == (xi2 < xii2 ? xi2 : xii2);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) && x3[1] == (xi2 > xii2 ? xi2 : xii2);

    x = x2.length();
    success += abs(x - sqrt(xii1*xii1 + xii2*xii2)) < 1.0e-15;

    x = length(x2);
    success += abs(x - sqrt(xii1*xii1 + xii2*xii2)) < 1.0e-15;

    while(sqr_length(x1) == 0.0)
    {
      xi1 = my_random();
      xi2 = my_random();      
      x1 = Vec2d(xi1, xi2);
    }
    x3 = normalize(x1);
    success += abs(x3[0] - xi1/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15 
               && abs(x3[1] - xi2/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15;

    x1.normalize();
    success += abs(x1[0] - xi1/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15 
               && abs(x1[1] - xi2/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15;
    
    x3 = orthogonal(x2);
    success += dot(x2, x3) == 0;

    x1 = Vec2d(xi1, xi2);
    x = cross(x1, x2);
    success += x == xi1*xii2 - xi2*xii1;

    while(x == 0)
    {
      xi1 = my_random();
      xi2 = my_random();
      x1 = Vec2d(xi1, xi2);
      x = cross(x1, x2);      
    }
    double xiii1 = my_random(), xiii2 = my_random();
    x3 = Vec2d(xiii1, xiii2);
    double y, z;
    linear_combine(x1, x2, x3, y, z);
    success += abs(y - (xii2*xiii1 - xii1*xiii2)/x) < 1.0e-15
               && abs(z - (xi1*xiii2 - xi2*xiii1)/x) < 1.0e-15;
  }
  if(success != 11)
  {
    cout << "Failure in test of Vec2d Vector operations";
    system ("pause");
    return 1;
  }
  success = 0;
  
  ////////////////////////////////////////////////
  //                V e c 2 f
  ////////////////////////////////////////////////

  // Constructors

  {
    Vec2f x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]);

    float xi1 = my_random(), xi2 = my_random();
    Vec2f x2(xi1, xi2);
    success += x2[0] == xi1 && x2[1] == xi2;

    int xii1 = gel_rand(), xii2 = gel_rand();
    Vec2i i1(xii1, xii2);
    Vec2f x3(i1);
    success += x3[0] == xii1 && x3[1] == xii2;
      
 // Doubles and floats cannot represent exactly the same numbers.
//    double xiii1 = my_random(), xiii2 = my_random();
//    Vec2d d1(xiii1, xiii2);
//    Vec2f x4(d1);
//    success += x4[0] == xiii1 && x4[1] == xiii2;

    float xiiii = my_random();
    Vec2f x5(xiiii);
    success += x5[0] == xiiii && x5[1] == xiiii;
  }
  if(success != 4)
  {
    cout << "Failure in test of Vec2f Constructors";
    system ("pause");
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec2f x1;    

    success += x1.get_dim() == 2;

    float xi1 = my_random(), xi2 = my_random();
    x1.set(xi1, xi2);
    success += x1[0] == xi1 && x1[1] == xi2;

    success += x1.get() == &x1[0];

    float temp = BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	temp = x1[i];
    success += temp == x1.min_coord();

    temp = -BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec2f Data manipulation";
    system ("pause");
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    float xi1 = my_random(), xi2 = my_random();
    while(xi1 == xi2)
    {
      xi1 = my_random();
      xi2 = my_random();
    }

    Vec2f x1(xi1, xi2), x2(xi1, xi2), x3(xi1), x4(xi2);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);

    if(xi1 < xi2)
    {
      success += x3.all_l(x4);
      success += !(x3.all_l(x2));
      success += x3.all_le(x2);
      success += !(x4.all_le(x3));
      success += x4.all_g(x3);
      success += !(x4.all_g(x2));
      success += x4.all_ge(x2);
      success += !(x3.all_ge(x4));
    }
    else 
    {
      success += x4.all_l(x3);
      success += !(x4.all_l(x2));
      success += x4.all_le(x2);
      success += !(x3.all_le(x4));
      success += x3.all_g(x4);
      success += !(x3.all_g(x2));
      success += x3.all_ge(x2);
      success += !(x4.all_ge(x3));
    }
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec2f Comparison operators";
    system ("pause");
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    float xi1 = my_random(), xi2 = my_random();
    Vec2f x1(xi1, xi2);
    float xi3 = my_random();
    x1 *= xi3;
    success += abs(x1[0] - xi1*xi3) < 1.0e-15 && abs(x1[1] - xi2*xi3) < 1.0e-15;

    while(xi3 == 0)
      xi3 = my_random();
    x1 = Vec2f(xi1, xi2);
    x1 /= xi3;
    success += abs(x1[0] - xi1/xi3) < 1.0e-15 && abs(x1[1] - xi2/xi3) < 1.0e-15;

    x1 = Vec2f(xi1, xi2);
    x1 += xi3;
    success += x1[0] == xi1 + xi3 && x1[1] == xi2 + xi3;

    x1 = Vec2f(xi1, xi2);
    x1 -= xi3;
    success += x1[0] == xi1 - xi3 && x1[1] == xi2 - xi3;

    float xii1 = my_random(), xii2 = my_random();
    Vec2f x2(xii1, xii2);
    x1 = Vec2f(xi1, xi2);
    x2 *= x1;
    success += abs(x2[0] - xi1*xii1) < 1.0e-15 && abs(x2[1] - xi2*xii2) < 1.0e-15;

    while(xi1 == 0)
      xi1 = my_random();
    while(xi2 == 0)
      xi2 = my_random();
    x1 = Vec2f(xi1, xi2);
    x2 = Vec2f(xii1, xii2);
    x2 /= x1;
    success += abs(x2[0] - xii1/xi1) < 1.0e-15 && abs(x2[1] - xii2/xi2) < 1.0e-15;

    x2 = Vec2f(xii1, xii2);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2;
    
    x2 = Vec2f(xii1, xii2);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec2f Assignment operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Unary operators

  {
    float xi1 = my_random(), xi2 = my_random();
    Vec2f x1 = -Vec2f(xi1, xi2);    
    success += x1[0] == -xi1 && x1[1] == -xi2;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec2f Unary operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Binary operators

  {
    float xi1 = my_random(), xi2 = my_random();
    Vec2f x1(xi1, xi2);
    float xii1 = my_random(), xii2 = my_random();
    while(xii1 == 0)
      xii1 = my_random();
    while(xii2 == 0)
      xii2 = my_random();
    Vec2f x2(xii1, xii2);
    Vec2f x3 = x1*x2;
    success += abs(x3[0] - xi1*xii1) < 1.0e-15 && abs(x3[1] - xi2*xii2) < 1.0e-15;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2;
    
    x3 = x1/x2;
    success += abs(x3[0] - xi1/xii1) < 1.0e-15 && abs(x3[1] - xi2/xii2) < 1.0e-15;

    float xi3 = my_random();
    x3 = x1*xi3;
    success += abs(x3[0] - xi1*xi3) < 1.0e-15 && abs(x3[1] - xi2*xi3) < 1.0e-15;
    
    x3 = xi3*x1;
    success += abs(x3[0] - xi1*xi3) < 1.0e-15 && abs(x3[1] - xi2*xi3) < 1.0e-15;

    float xi4 = my_random();
    x3 = xi4*x1;
    success += abs(x3[0] - xi1*xi4) < 1-0e-15 && abs(x3[1] - xi2*xi4) < 1.0e-15;
    
    x3 = x1*xi4;
    success += abs(x3[0] - xi1*xi4) < 1-0e-15 && abs(x3[1] - xi2*xi4) < 1.0e-15;

    int xi5 = gel_rand();
    x3 = xi5*x1;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 && abs(x3[1] - xi2*xi5) < 1.0e-15;

    x3 = x1*xi5;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 && abs(x3[1] - xi2*xi5) < 1.0e-15;
    
    while(xi3 == 0)
      xi3 = my_random();
    x3 = x1/xi3;
    success += abs(x3[0] - xi1/xi3) < 1.0e-15 && abs(x3[1] - xi2/xi3) < 1.0e-15;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec2f Binary operators";
    system ("pause");
    return 1;
  }
  success = 0;

  // Vector operations

  {
    float xi1 = my_random(), xi2 = my_random();
    Vec2f x1(xi1, xi2);
    float xii1 = my_random(), xii2 = my_random();
    Vec2f x2(xii1, xii2);
    float x = dot(x1, x2);
    success += (x - xi1*xii1 - xi2*xii2) < 1.0e-7;

    x = sqr_length(x1);
    success += abs(x - xi1*xi1 - xi2*xi2) < 1.0e-15;

    Vec2f x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) && x3[1] == (xi2 < xii2 ? xi2 : xii2);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) && x3[1] == (xi2 > xii2 ? xi2 : xii2);

    x = x2.length();
    success += abs(x - sqrt(xii1*xii1 + xii2*xii2)) < 1.0e-15;

    x = length(x2);
    success += abs(x - sqrt(xii1*xii1 + xii2*xii2)) < 1.0e-15;

    while(sqr_length(x1) == 0.0)
    {
      xi1 = my_random();
      xi2 = my_random();      
      x1 = Vec2f(xi1, xi2);
    }
    x3 = normalize(x1);
    success += abs(x3[0] - xi1/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15 
               && abs(x3[1] - xi2/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15;

    x1.normalize();
    success += abs(x1[0] - xi1/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15 
               && abs(x1[1] - xi2/sqrt(xi1*xi1 + xi2*xi2)) < 1.0e-15;
    
    x3 = orthogonal(x2);
    success += abs(dot(x2, x3)) < 1.0e-15;

    x1 = Vec2f(xi1, xi2);
    x = cross(x1, x2);
    success += abs(x - xi1*xii2 + xi2*xii1) < 1.0e-15;

    while(x == 0)
    {
      xi1 = my_random();
      xi2 = my_random();
      x1 = Vec2f(xi1, xi2);
      x = cross(x1, x2);      
    }
    float xiii1 = my_random(), xiii2 = my_random();
    x3 = Vec2f(xiii1, xiii2);
    float y, z;
    linear_combine(x1, x2, x3, y, z);
    success += abs(y - (xii2*xiii1 - xii1*xiii2)/x) < 1.0e-15
               && abs(z - (xi1*xiii2 - xi2*xiii1)/x) < 1.0e-15;
  }
  if(success != 11)
  {
    cout << "Failure in test of Vec2f Vector operations";
    system ("pause");
    return 1;
  }
  success = 0;

  ////////////////////////////////////////////////
  //                V e c 2 i
  ////////////////////////////////////////////////

  // Constructors

  {
/* Vec2i default initialization ?
 *
    Vec2i x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]);
*/

    int xi1 = gel_rand(), xi2 = gel_rand();
    Vec2i x2(xi1, xi2);
    success += x2[0] == xi1 && x2[1] == xi2;

/* Constuctor non-existent !
 *
    int xii1 = my_random(), xii2 = my_random();
    Vec2d d1(xii1, xii2);
    Vec2i x3(i1);
    success += x3[0] == xii1 && x3[1] == xii2;
*/

    float xiii1 = my_random(), xiii2 = my_random();
    Vec2f f1(xiii1, xiii2);
    Vec2i x4(f1);
    success += x4[0] == static_cast<int>(xiii1) && x4[1] == static_cast<int>(xiii2);

/* Constuctor non-existent !
 *
    int xiiii = gel_rand();
    Vec2i x5(xiiii);
    success += x5[0] == xiiii && x5[1] == xiiii;
*/
  }
  if(success != 2)
  {
    cout << "Failure in test of Vec2i Constructors";
    system ("pause");
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec2i x1;    

    success += x1.get_dim() == 2;

    int xi1 = gel_rand(), xi2 = gel_rand();
    x1.set(xi1, xi2);
    success += x1[0] == xi1 && x1[1] == xi2;

    success += x1.get() == &x1[0];

    int temp = INT_MAX;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	temp = x1[i];
    success += temp == x1.min_coord();

    temp = -INT_MAX;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec2i Data manipulation";
    system ("pause");
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand();
    while(xi1 == xi2)
    {
      xi1 = gel_rand();
      xi2 = gel_rand();
    }

    Vec2i x1(xi1, xi2), x2(xi1, xi2), x3(xi1, xi1), x4(xi2, xi2);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);

    if(xi1 < xi2)
    {
      success += x3.all_l(x4);
      success += !(x3.all_l(x2));
      success += x3.all_le(x2);
      success += !(x4.all_le(x3));
      success += x4.all_g(x3);
      success += !(x4.all_g(x2));
      success += x4.all_ge(x2);
      success += !(x3.all_ge(x4));
    }
    else 
    {
      success += x4.all_l(x3);
      success += !(x4.all_l(x2));
      success += x4.all_le(x2);
      success += !(x3.all_le(x4));
      success += x3.all_g(x4);
      success += !(x3.all_g(x2));
      success += x3.all_ge(x2);
      success += !(x4.all_ge(x3));
    }
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec2i Comparison operators";
    system ("pause");
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand();
    Vec2i x1(xi1, xi2);
    int xi3 = gel_rand();
    x1 *= xi3;
    success += x1[0] == xi1*xi3 && x1[1] == xi2*xi3;

    while(xi3 == 0)
      xi3 = gel_rand();
    x1 = Vec2i(xi1, xi2);
    x1 /= xi3;
    success += x1[0] == xi1/xi3 && x1[1] == xi2/xi3;

    x1 = Vec2i(xi1, xi2);
    x1 += xi3;
    success += x1[0] == xi1 + xi3 && x1[1] == xi2 + xi3;

    x1 = Vec2i(xi1, xi2);
    x1 -= xi3;
    success += x1[0] == xi1 - xi3 && x1[1] == xi2 - xi3;

    int xii1 = gel_rand(), xii2 = gel_rand();
    Vec2i x2(xii1, xii2);
    x1 = Vec2i(xi1, xi2);
    x2 *= x1;
    success += x2[0] == xi1*xii1 && x2[1] == xi2*xii2;

    while(xi1 == 0)
      xi1 = gel_rand();
    while(xi2 == 0)
      xi2 = gel_rand();
    x1 = Vec2i(xi1, xi2);
    x2 = Vec2i(xii1, xii2);
    x2 /= x1;
    success += x2[0] == xii1/xi1 && x2[1] == xii2/xi2;

    x2 = Vec2i(xii1, xii2);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2;
    
    x2 = Vec2i(xii1, xii2);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec2i Assignment operators";
    return 1;
  }
  success = 0;

  // Unary operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand();
    Vec2i x1 = -Vec2i(xi1, xi2);    
    success += x1[0] == -xi1 && x1[1] == -xi2;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec2i Unary operators";
    return 1;
  }
  success = 0;

  // Binary operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand();
    Vec2i x1(xi1, xi2);
    int xii1 = gel_rand(), xii2 = gel_rand();
    while(xii1 == 0)
      xii1 = gel_rand();
    while(xii2 == 0)
      xii2 = gel_rand();
    Vec2i x2(xii1, xii2);
    Vec2i x3 = x1*x2;
    success += x3[0] == xi1*xii1 && x3[1] == xi2*xii2;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2;
    
    x3 = x1/x2;
    success += x3[0] == xi1/xii1 && x3[1] == xi2/xii2;

    int xi3 = gel_rand();
    x3 = x1*xi3;
    success += x3[0] == xi1*xi3 && x3[1] == xi2*xi3;
    
    x3 = xi3*x1;
    success += x3[0] == xi1*xi3 && x3[1] == xi2*xi3;

    float xi4 = INT_MAX*my_random();
    x3 = xi4*x1;
    success += x3[0] == xi1*static_cast<int>(xi4) && x3[1] == xi2*static_cast<int>(xi4);
    
    x3 = x1*xi4;
    success += x3[0] == xi1*static_cast<int>(xi4) && x3[1] == xi2*static_cast<int>(xi4);

    double xi5 = INT_MAX*my_random();
    x3 = xi5*x1;
    success += x3[0] == xi1*static_cast<int>(xi5) && x3[1] == xi2*static_cast<int>(xi5);
    
    x3 = x1*xi5;
    success += x3[0] == xi1*static_cast<int>(xi5) && x3[1] == xi2*static_cast<int>(xi5);

    while(xi3 == 0)
      xi3 = gel_rand();
    x3 = x1/xi3;
    success += x3[0] == xi1/xi3 && x3[1] == xi2/xi3;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec2i Binary operators";
    return 1;
  }
  success = 0;

  // Vector operations

  {
    int xi1 = gel_rand(), xi2 = gel_rand();
    Vec2i x1(xi1, xi2);
    int xii1 = gel_rand(), xii2 = gel_rand();
    Vec2i x2(xii1, xii2);
    int x = dot(x1, x2);
    success += x == xi1*xii1 + xi2*xii2;

    x = sqr_length(x1);
    success += x == xi1*xi1 + xi2*xi2;

    Vec2i x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) && x3[1] == (xi2 < xii2 ? xi2 : xii2);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) && x3[1] == (xi2 > xii2 ? xi2 : xii2);
  }
  if(success != 4)
  {
    cout << "Failure in test of Vec2i Vector operations";
    return 1;
  }
  success = 0;

  ////////////////////////////////////////////////
  //                V e c 3 i
  ////////////////////////////////////////////////

  // Constructors

  {
/* Vec3i default initialization ?
 *
    Vec3i x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]) && CGLA::isnan(x1[2]);
*/

    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    Vec3i x2(xi1, xi2, xi3);
    success += x2[0] == xi1 && x2[1] == xi2 && x2[2] == xi3;

    int xiiii = gel_rand();
    Vec3i x3(xiiii);
    success += x3[0] == xiiii && x3[1] == xiiii && x3[2] == xiiii;

/* Constuctor non-existent !
 *
    int xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3d d1(xii1, xii2, xii3);
    Vec3i x4(d1);
    success += x4[0] == xii1 && x4[1] == xii2 && x4[2] == xii3;
*/

    float xiii1 = my_random(), xiii2 = my_random(), xiii3 = my_random();
    Vec3f f1(xiii1, xiii2, xiii3);
    Vec3i x5(f1);
    success += x5[0] == static_cast<int>(xiii1) 
               && x5[1] == static_cast<int>(xiii2)
               && x5[2] == static_cast<int>(xiii3);

    unsigned char xiiii1 = 256*my_random(), 
                  xiiii2 = 256*my_random(), 
                  xiiii3 = 256*my_random();
    Vec3uc uc1(xiiii1, xiiii2, xiiii3);
    Vec3i x6(uc1);
    success += x6[0] == xiiii1 && x6[1] == xiiii2 && x6[2] == xiiii3;

    unsigned short int xiiv1 = USHRT_MAX*my_random(),
                       xiiv2 = USHRT_MAX*my_random(),
                       xiiv3 = USHRT_MAX*my_random();
    Vec3usi usi1(xiiv1, xiiv2, xiiv3);
    Vec3i x7(usi1);
    success += x7[0] == xiiv1 && x7[1] == xiiv2 && x7[2] == xiiv3;
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3i Constructors";
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec3i x1;    

    success += x1.get_dim() == 3;

    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    x1.set(xi1, xi2, xi3);
    success += x1[0] == xi1 && x1[1] == xi2 && x1[2] == xi3;

    success += x1.get() == &x1[0];

    int temp = INT_MAX;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	temp = x1[i];
    success += temp == x1.min_coord();

    temp = -INT_MAX;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3i Data manipulation";
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    while(xi1 >= xi2)
    {
      xi1 = gel_rand();
      xi2 = gel_rand();
    }
    while(xi3 <= xi2)
      xi3 = gel_rand();

    Vec3i x1(xi1, xi2, xi3), x2(xi1, xi2, xi3), 
          x3(xi1, xi1, xi1), x4(xi2, xi2, xi2), x5(xi3, xi3, xi3);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);
    success += x3.all_l(x4);
    success += !(x3.all_l(x2));
    success += x3.all_le(x2);
    success += !(x4.all_le(x2));
    success += x4.all_g(x3);
    success += !(x4.all_g(x2));
    success += x5.all_ge(x2);
    success += !(x4.all_ge(x2));
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec3i Comparison operators";
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    Vec3i x1(xi1, xi2, xi3);
    int xii = gel_rand();
    x1 *= xii;
    success += x1[0] == xi1*xii && x1[1] == xi2*xii && x1[2] == xi3*xii;

    while(xii == 0)
      xii = gel_rand();
    x1 = Vec3i(xi1, xi2, xi3);
    x1 /= xii;
    success += x1[0] == xi1/xii && x1[1] == xi2/xii && x1[2] == xi3/xii;

    x1 = Vec3i(xi1, xi2, xi3);
    x1 += xii;
    success += x1[0] == xi1 + xii && x1[1] == xi2 + xii && x1[2] == xi3 + xii;

    x1 = Vec3i(xi1, xi2, xi3);
    x1 -= xii;
    success += x1[0] == xi1 - xii && x1[1] == xi2 - xii && x1[2] == xi3 - xii;

    int xii1 = gel_rand(), xii2 = gel_rand(), xii3 = gel_rand();
    Vec3i x2(xii1, xii2, xii3);
    x1 = Vec3i(xi1, xi2, xi3);
    x2 *= x1;
    success += x2[0] == xi1*xii1 && x2[1] == xi2*xii2 && x2[2] == xi3*xii3;

    while(xi1 == 0)
      xi1 = gel_rand();
    while(xi2 == 0)
      xi2 = gel_rand();
    while(xi3 == 0)
      xi3 = gel_rand();
    x1 = Vec3i(xi1, xi2, xi3);
    x2 = Vec3i(xii1, xii2, xii3);
    x2 /= x1;
    success += x2[0] == xii1/xi1 && x2[1] == xii2/xi2 && x2[2] == xii3/xi3;

    x2 = Vec3i(xii1, xii2, xii3);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2 && x2[2] == xii3 + xi3;
    
    x2 = Vec3i(xii1, xii2, xii3);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2 && x2[2] == xii3 - xi3;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec3i Assignment operators";
    return 1;
  }
  success = 0;

  // Unary operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    Vec3i x1 = -Vec3i(xi1, xi2, xi3);    
    success += x1[0] == -xi1 && x1[1] == -xi2 && x1[2] == -xi3;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec3i Unary operators";
    return 1;
  }
  success = 0;

  // Binary operators

  {
    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    Vec3i x1(xi1, xi2, xi3);
    int xii1 = gel_rand(), xii2 = gel_rand(), xii3 = gel_rand();
    while(xii1 == 0)
      xii1 = gel_rand();
    while(xii2 == 0)
      xii2 = gel_rand();
    while(xii3 == 0)
      xii3 = gel_rand();
    Vec3i x2(xii1, xii2, xii3);
    Vec3i x3 = x1*x2;
    success += x3[0] == xi1*xii1 && x3[1] == xi2*xii2 && x3[2] == xi3*xii3;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2 && x3[2] == xi3 + xii3;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2 && x3[2] == xi3 - xii3;
    
    x3 = x1/x2;
    success += x3[0] == xi1/xii1 && x3[1] == xi2/xii2 && x3[2] == xi3/xii3;

    int xii = gel_rand();
    x3 = x1*xii;
    success += x3[0] == xi1*xii && x3[1] == xi2*xii && x3[2] == xi3*xii;
    
    x3 = xii*x1;
    success += x3[0] == xi1*xii && x3[1] == xi2*xii && x3[2] == xi3*xii;

    float xi4 = INT_MAX*my_random();
    x3 = xi4*x1;
    success += x3[0] == xi1*static_cast<int>(xi4) 
               && x3[1] == xi2*static_cast<int>(xi4) 
               && x3[2] == xi3*static_cast<int>(xi4);
    
    x3 = x1*xi4;
    success += x3[0] == xi1*static_cast<int>(xi4) 
               && x3[1] == xi2*static_cast<int>(xi4) 
               && x3[2] == xi3*static_cast<int>(xi4);

    double xi5 = INT_MAX*my_random();
    x3 = xi5*x1;
    success += x3[0] == xi1*static_cast<int>(xi5) 
               && x3[1] == xi2*static_cast<int>(xi5) 
               && x3[2] == xi3*static_cast<int>(xi5);
    
    x3 = x1*xi5;
    success += x3[0] == xi1*static_cast<int>(xi5) 
               && x3[1] == xi2*static_cast<int>(xi5) 
               && x3[2] == xi3*static_cast<int>(xi5);

    while(xii == 0)
      xii = gel_rand();
    x3 = x1/xii;
    success += x3[0] == xi1/xii && x3[1] == xi2/xii && x3[2] == xi3/xii;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec3i Binary operators";
    return 1;
  }
  success = 0;

  // Vector operations

  {
    int xi1 = gel_rand(), xi2 = gel_rand(), xi3 = gel_rand();
    Vec3i x1(xi1, xi2, xi3);
    int xii1 = gel_rand(), xii2 = gel_rand(), xii3 = gel_rand();
    Vec3i x2(xii1, xii2, xii3);
    int x = dot(x1, x2);
    success += x == xi1*xii1 + xi2*xii2 + xi3*xii3;

    x = sqr_length(x1);
    success += x == xi1*xi1 + xi2*xi2 + xi3*xi3;

    Vec3i x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 < xii2 ? xi2 : xii2)
               && x3[2] == (xi3 < xii3 ? xi3 : xii3);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 > xii2 ? xi2 : xii2)
               && x3[2] == (xi3 > xii3 ? xi3 : xii3);

    x3 = cross(x1, x2);
    success += x3[0] == xi2*xii3 - xi3*xii2
               && x3[1] == xi3*xii1 - xi1*xii3
               && x3[2] == xi1*xii2 - xi2*xii1;
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3i Vector operations";
    return 1;
  }
  success = 0;

  ////////////////////////////////////////////////
  //                V e c 3 f
  ////////////////////////////////////////////////

  // Constructors

  {
    Vec3f x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]) && CGLA::isnan(x1[2]);

    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3f x2(xi1, xi2, xi3);
    success += x2[0] == xi1 && x2[1] == xi2 && x2[2] == xi3;

    float xiiii = my_random();
    Vec3f x3(xiiii);
    success += x3[0] == xiiii && x3[1] == xiiii && x3[2] == xiiii;

      // Doubles and floats cannot represent exactly the same numbers.
//    double xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
//    Vec3d d1(xii1, xii2, xii3);
//    Vec3f x4(d1);
//    success += x4[0] == xii1 && x4[1] == xii2 && x4[2] == xii3;

    int xiii1 = gel_rand(), xiii2 = gel_rand(), xiii3 = gel_rand();
    Vec3i i1(xiii1, xiii2, xiii3);
    Vec3f x5(i1);
    success += x5[0] == xiii1 && x5[1] == xiii2 && x5[2] == xiii3;

    float xiiii1 = my_random(), xiiii2 = my_random(), xiiii3 = my_random();
    Quatf q(xiiii1, xiiii2, xiiii3, 1.0);
    Vec3f x6 = q.get_imaginary_part();
    success += x6[0] == xiiii1 && x6[1] == xiiii2 && x6[2] == xiiii3;

    unsigned short int xiiv1 = USHRT_MAX*my_random(),
                       xiiv2 = USHRT_MAX*my_random(),
                       xiiv3 = USHRT_MAX*my_random();
    Vec3usi usi1(xiiv1, xiiv2, xiiv3);
    Vec3f x7(usi1);
    success += x7[0] == xiiv1 && x7[1] == xiiv2 && x7[2] == xiiv3;
  }
  if(success != 7)
  {
    cout << "Failure in test of Vec3f Constructors";
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec3f x1;    

    success += x1.get_dim() == 3;

    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    x1.set(xi1, xi2, xi3);
    success += x1[0] == xi1 && x1[1] == xi2 && x1[2] == xi3;

    success += x1.get() == &x1[0];

    float temp = BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	temp = x1[i];
    success += temp == x1.min_coord();

    temp = -BIG;
    for(unsigned int i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3f Data manipulation";
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    while(xi1 >= xi2)
    {
      xi1 = my_random();
      xi2 = my_random();
    }
    while(xi3 <= xi2)
      xi3 = my_random();

    Vec3f x1(xi1, xi2, xi3), x2(xi1, xi2, xi3), 
          x3(xi1, xi1, xi1), x4(xi2, xi2, xi2), x5(xi3, xi3, xi3);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);
    success += x3.all_l(x4);
    success += !(x3.all_l(x2));
    success += x3.all_le(x2);
    success += !(x4.all_le(x2));
    success += x4.all_g(x3);
    success += !(x4.all_g(x2));
    success += x5.all_ge(x2);
    success += !(x4.all_ge(x2));
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec3f Comparison operators";
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3f x1(xi1, xi2, xi3);
    float xii = my_random();
    x1 *= xii;
    success += abs(x1[0] - xi1*xii) < 1.0e-15 
               && abs(x1[1] - xi2*xii) < 1.0e-15 
               && abs(x1[2] - xi3*xii) < 1.0e-15;

    while(xii == 0)
      xii = my_random();
    x1 = Vec3f(xi1, xi2, xi3);
    x1 /= xii;
    success += abs(x1[0] - xi1/xii) < 1.0e-15 
               && abs(x1[1] - xi2/xii) < 1.0e-15 
               && abs(x1[2] - xi3/xii) < 1.0e-15;

    x1 = Vec3f(xi1, xi2, xi3);
    x1 += xii;
    success += x1[0] == xi1 + xii && x1[1] == xi2 + xii && x1[2] == xi3 + xii;

    x1 = Vec3f(xi1, xi2, xi3);
    x1 -= xii;
    success += x1[0] == xi1 - xii && x1[1] == xi2 - xii && x1[2] == xi3 - xii;

    float xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3f x2(xii1, xii2, xii3);
    x1 = Vec3f(xi1, xi2, xi3);
    x2 *= x1;
    success += abs(x2[0] - xi1*xii1) < 1.0e-15 
               && abs(x2[1] - xi2*xii2) < 1.0e-15 
               && abs(x2[2] - xi3*xii3) < 1.0e-15;

    while(xi1 == 0)
      xi1 = my_random();
    while(xi2 == 0)
      xi2 = my_random();
    while(xi3 == 0)
      xi3 = my_random();
    x1 = Vec3f(xi1, xi2, xi3);
    x2 = Vec3f(xii1, xii2, xii3);
    x2 /= x1;
    success += abs(x2[0] - xii1/xi1) < 1.0e-15 
               && abs(x2[1] - xii2/xi2) < 1.0e-15 
               && abs(x2[2] - xii3/xi3) < 1.0e-15;

    x2 = Vec3f(xii1, xii2, xii3);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2 && x2[2] == xii3 + xi3;
    
    x2 = Vec3f(xii1, xii2, xii3);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2 && x2[2] == xii3 - xi3;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec3f Assignment operators";
    return 1;
  }
  success = 0;

  // Unary operators

  {
    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3f x1 = -Vec3f(xi1, xi2, xi3);    
    success += x1[0] == -xi1 && x1[1] == -xi2 && x1[2] == -xi3;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec3f Unary operators";
    return 1;
  }
  success = 0;

  // Binary operators

  {
    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3f x1(xi1, xi2, xi3);
    float xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    while(xii1 == 0)
      xii1 = my_random();
    while(xii2 == 0)
      xii2 = my_random();
    while(xii3 == 0)
      xii3 = my_random();
    Vec3f x2(xii1, xii2, xii3);
    Vec3f x3 = x1*x2;
    success += abs(x3[0] - xi1*xii1) < 1.0e-15 
               && abs(x3[1] - xi2*xii2) < 1.0e-15 
               && abs(x3[2] - xi3*xii3) < 1.0e-15;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2 && x3[2] == xi3 + xii3;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2 && x3[2] == xi3 - xii3;
    
    x3 = x1/x2;
    success += abs(x3[0] - xi1/xii1) < 1.0e-15 
               && abs(x3[1] - xi2/xii2) < 1.0e-15 
               && abs(x3[2] - xi3/xii3) < 1.0e-15;

    float xii = my_random();
    x3 = x1*xii;
    success += abs(x3[0] - xi1*xii) < 1.0e-15
               && abs(x3[1] - xi2*xii) < 1.0e-15 
               && abs(x3[2] - xi3*xii) < 1.0e-15;
    
    x3 = xii*x1;
    success += abs(x3[0] - xi1*xii) < 1.0e-15 
               && abs(x3[1] - xi2*xii) < 1.0e-15 
               && abs(x3[2] - xi3*xii) < 1.0e-15;

    int xi4 = gel_rand();
    x3 = xi4*x1;
    success += abs(x3[0] - xi1*xi4) < 1.0e-15 
               && abs(x3[1] - xi2*xi4) < 1.0e-15 
               && abs(x3[2] - xi3*xi4) < 1.0e-15;
    
    x3 = x1*xi4;
    success += abs(x3[0] - xi1*xi4) < 1.0e-15 
               && abs(x3[1] - xi2*xi4) < 1.0e-15 
               && abs(x3[2] - xi3*xi4) < 1.0e-15;

    double xi5 = my_random();
    x3 = xi5*x1;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 
               && abs(x3[1] - xi2*xi5) < 1.0e-15 
               && abs(x3[2] - xi3*xi5) < 1.0e-15;
    
    x3 = x1*xi5;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 
               && abs(x3[1] - xi2*xi5) < 1.0e-15 
               && abs(x3[2] - xi3*xi5) < 1.0e-15;

    while(xii == 0)
      xii = my_random();
    x3 = x1/xii;
    success += abs(x3[0] - xi1/xii) < 1.0e-15 
               && abs(x3[1] - xi2/xii) < 1.0e-15 
               && abs(x3[2] - xi3/xii) < 1.0e-15;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec3f Binary operators";
    return 1;
  }
  success = 0;

  // Vector operations

  {
    float xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3f x1(xi1, xi2, xi3);
    float xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3f x2(xii1, xii2, xii3);
    float x = dot(x1, x2);
    success += abs(x - xi1*xii1 - xi2*xii2 - xi3*xii3) < 1.0e-15;

    x = sqr_length(x1);
    success += abs(x - xi1*xi1 - xi2*xi2 - xi3*xi3) < 1.0e-15;

    Vec3f x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 < xii2 ? xi2 : xii2)
               && x3[2] == (xi3 < xii3 ? xi3 : xii3);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 > xii2 ? xi2 : xii2)
               && x3[2] == (xi3 > xii3 ? xi3 : xii3);

    Vec3f x4;
    orthogonal(x2, x3, x4);
    success += abs(dot(x2, x3)) < 1.0e-15;

    x3 = cross(x1, x2);
    success += abs(x3[0] - xi2*xii3 + xi3*xii2) < 1.0e-15
               && abs(x3[1] - xi3*xii1 + xi1*xii3) < 1.0e-15
               && abs(x3[2] - xi1*xii2 + xi2*xii1) < 1.0e-15;

    float theta, phi, r;
    x1.get_spherical(theta, phi, r);
    success += abs(theta - acos(xi3/sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3))) < 1.0e-15
               && abs(phi - atan(xi2/xi1)) < 1.0e-15 
               && abs(r - sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3)) < 1.0e-15;

    x3.set_spherical(theta, phi, r);
    success += abs(x3[0] - xi1) < 1.0e-15 
               && abs(x3[1] - xi2) < 1.0e-15 
               && abs(x3[2] - xi3) < 1.0e-15;
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec3f Vector operations";
    return 1;
  }
  success = 0;

  ////////////////////////////////////////////////
  //                V e c 3 d
  ////////////////////////////////////////////////

  // Constructors

  {
    Vec3d x1;
    success += CGLA::isnan(x1[0]) && CGLA::isnan(x1[1]) && CGLA::isnan(x1[2]);

    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3d x2(xi1, xi2, xi3);
    success += x2[0] == xi1 && x2[1] == xi2 && x2[2] == xi3;

    double xiiii = my_random();
    Vec3d x3(xiiii);
    success += x3[0] == xiiii && x3[1] == xiiii && x3[2] == xiiii;

    float xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3f f1(xii1, xii2, xii3);
    Vec3d x4(f1);
    success += x4[0] == xii1 && x4[1] == xii2 && x4[2] == xii3;

    int xiii1 = gel_rand(), xiii2 = gel_rand(), xiii3 = gel_rand();
    Vec3i i1(xiii1, xiii2, xiii3);
    Vec3d x5(i1);
    success += x5[0] == xiii1 && x5[1] == xiii2 && x5[2] == xiii3;

 
    double xiiii1 = my_random(), xiiii2 = my_random(), xiiii3 = my_random();
    Quatd q(xiiii1, xiiii2, xiiii3, 1.0);
    Vec3d x6 = q.get_imaginary_part();
    success += x6[0] == xiiii1 && x6[1] == xiiii2 && x6[2] == xiiii3;

    unsigned short int xiiv1 = USHRT_MAX*my_random(),
                       xiiv2 = USHRT_MAX*my_random(),
                       xiiv3 = USHRT_MAX*my_random();
    Vec3usi usi1(xiiv1, xiiv2, xiiv3);
    Vec3d x7(usi1);
    success += x7[0] == xiiv1 && x7[1] == xiiv2 && x7[2] == xiiv3;

  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3d Constructors";
    return 1;
  }
  success = 0;

  // Data manipulation

  {
    Vec3d x1;    

    success += x1.get_dim() == 3;

    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    x1.set(xi1, xi2, xi3);
    success += x1[0] == xi1 && x1[1] == xi2 && x1[2] == xi3;

    success += x1.get() == &x1[0];

    double temp = BIG;
    for(double i = 0; i < x1.get_dim(); ++i)
      if(x1[i] < temp)
	temp = x1[i];
    success += temp == x1.min_coord();

    temp = -BIG;
    for(double i = 0; i < x1.get_dim(); ++i)
      if(x1[i] > temp)
	temp = x1[i];
    success += temp == x1.max_coord();
  }
  if(success != 5)
  {
    cout << "Failure in test of Vec3d Data manipulation";
    return 1;
  }
  success = 0;

  // Comparison operators

  {
    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    while(xi1 >= xi2)
    {
      xi1 = my_random();
      xi2 = my_random();
    }
    while(xi3 <= xi2)
      xi3 = my_random();

    Vec3d x1(xi1, xi2, xi3), x2(xi1, xi2, xi3), 
          x3(xi1, xi1, xi1), x4(xi2, xi2, xi2), x5(xi3, xi3, xi3);    
    success += x1 == x2;
    success += !(x1 == x3);
    success += x3 == xi1;
    success += !(x3 == xi2);
    success += x1 != x3;
    success += !(x1 != x2);
    success += x3.all_l(x4);
    success += !(x3.all_l(x2));
    success += x3.all_le(x2);
    success += !(x4.all_le(x2));
    success += x4.all_g(x3);
    success += !(x4.all_g(x2));
    success += x5.all_ge(x2);
    success += !(x4.all_ge(x2));
  }
  if(success != 14)
  {
    cout << "Failure in test of Vec3d Comparison operators";
    return 1;
  }
  success = 0;
  
  // Assignment operators

  {
    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3d x1(xi1, xi2, xi3);
    double xii = my_random();
    x1 *= xii;
    success += abs(x1[0] - xi1*xii) < 1.0e-15 
               && abs(x1[1] - xi2*xii) < 1.0e-15 
               && abs(x1[2] - xi3*xii) < 1.0e-15;

    while(xii == 0)
      xii = my_random();
    x1 = Vec3d(xi1, xi2, xi3);
    x1 /= xii;
    success += abs(x1[0] - xi1/xii) < 1.0e-15 
               && abs(x1[1] - xi2/xii) < 1.0e-15 
               && abs(x1[2] - xi3/xii) < 1.0e-15;

    x1 = Vec3d(xi1, xi2, xi3);
    x1 += xii;
    success += x1[0] == xi1 + xii && x1[1] == xi2 + xii && x1[2] == xi3 + xii;

    x1 = Vec3d(xi1, xi2, xi3);
    x1 -= xii;
    success += x1[0] == xi1 - xii && x1[1] == xi2 - xii && x1[2] == xi3 - xii;

    double xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3d x2(xii1, xii2, xii3);
    x1 = Vec3d(xi1, xi2, xi3);
    x2 *= x1;
    success += abs(x2[0] - xi1*xii1) < 1.0e-15 
               && abs(x2[1] - xi2*xii2) < 1.0e-15 
               && abs(x2[2] - xi3*xii3) < 1.0e-15;

    while(xi1 == 0)
      xi1 = my_random();
    while(xi2 == 0)
      xi2 = my_random();
    while(xi3 == 0)
      xi3 = my_random();
    x1 = Vec3d(xi1, xi2, xi3);
    x2 = Vec3d(xii1, xii2, xii3);
    x2 /= x1;
    success += abs(x2[0] - xii1/xi1) < 1.0e-15 
               && abs(x2[1] - xii2/xi2) < 1.0e-15 
               && abs(x2[2] - xii3/xi3) < 1.0e-15;

    x2 = Vec3d(xii1, xii2, xii3);
    x2 += x1;
    success += x2[0] == xii1 + xi1 && x2[1] == xii2 + xi2 && x2[2] == xii3 + xi3;
    
    x2 = Vec3d(xii1, xii2, xii3);
    x2 -= x1;
    success += x2[0] == xii1 - xi1 && x2[1] == xii2 - xi2 && x2[2] == xii3 - xi3;    
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec3d Assignment operators";
    return 1;
  }
  success = 0;

  // Unary operators

  {
    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3d x1 = -Vec3d(xi1, xi2, xi3);    
    success += x1[0] == -xi1 && x1[1] == -xi2 && x1[2] == -xi3;
  }
  if(success != 1)
  {
    cout << "Failure in test of Vec3d Unary operators";
    return 1;
  }
  success = 0;

  // Binary operators

  {
    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3d x1(xi1, xi2, xi3);
    double xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    while(xii1 == 0)
      xii1 = my_random();
    while(xii2 == 0)
      xii2 = my_random();
    while(xii3 == 0)
      xii3 = my_random();
    Vec3d x2(xii1, xii2, xii3);
    Vec3d x3 = x1*x2;
    success += abs(x3[0] - xi1*xii1) < 1.0e-15 
               && abs(x3[1] - xi2*xii2) < 1.0e-15 
               && abs(x3[2] - xi3*xii3) < 1.0e-15;
    
    x3 = x1 + x2;
    success += x3[0] == xi1 + xii1 && x3[1] == xi2 + xii2 && x3[2] == xi3 + xii3;

    x3 = x1 - x2;
    success += x3[0] == xi1 - xii1 && x3[1] == xi2 - xii2 && x3[2] == xi3 - xii3;
    
    x3 = x1/x2;
    success += abs(x3[0] - xi1/xii1) < 1.0e-15 
               && abs(x3[1] - xi2/xii2) < 1.0e-15 
               && abs(x3[2] - xi3/xii3) < 1.0e-15;

    double xii = my_random();
    x3 = x1*xii;
    success += abs(x3[0] - xi1*xii) < 1.0e-15
               && abs(x3[1] - xi2*xii) < 1.0e-15 
               && abs(x3[2] - xi3*xii) < 1.0e-15;
    
    x3 = xii*x1;
    success += abs(x3[0] - xi1*xii) < 1.0e-15 
               && abs(x3[1] - xi2*xii) < 1.0e-15 
               && abs(x3[2] - xi3*xii) < 1.0e-15;

    int xi4 = gel_rand();
    x3 = xi4*x1;
    success += abs(x3[0] - xi1*xi4) < 1.0e-15 
               && abs(x3[1] - xi2*xi4) < 1.0e-15 
               && abs(x3[2] - xi3*xi4) < 1.0e-15;
    
    x3 = x1*xi4;
    success += abs(x3[0] - xi1*xi4) < 1.0e-15 
               && abs(x3[1] - xi2*xi4) < 1.0e-15 
               && abs(x3[2] - xi3*xi4) < 1.0e-15;

    float xi5 = my_random();
    x3 = xi5*x1;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 
               && abs(x3[1] - xi2*xi5) < 1.0e-15 
               && abs(x3[2] - xi3*xi5) < 1.0e-15;
    
    x3 = x1*xi5;
    success += abs(x3[0] - xi1*xi5) < 1.0e-15 
               && abs(x3[1] - xi2*xi5) < 1.0e-15 
               && abs(x3[2] - xi3*xi5) < 1.0e-15;

    while(xii == 0)
      xii = my_random();
    x3 = x1/xii;
    success += abs(x3[0] - xi1/xii) < 1.0e-15 
               && abs(x3[1] - xi2/xii) < 1.0e-15 
               && abs(x3[2] - xi3/xii) < 1.0e-15;
  } 
  if(success != 11)
  {
    cout << "Failure in test of Vec3d Binary operators";
    return 1;
  }
  success = 0;

  // Vector operations

  {
    double xi1 = my_random(), xi2 = my_random(), xi3 = my_random();
    Vec3d x1(xi1, xi2, xi3);
    double xii1 = my_random(), xii2 = my_random(), xii3 = my_random();
    Vec3d x2(xii1, xii2, xii3);
    double x = dot(x1, x2);
    success += abs(x - xi1*xii1 - xi2*xii2 - xi3*xii3) < 1.0e-15;

    x = sqr_length(x1);
    success += abs(x - xi1*xi1 - xi2*xi2 - xi3*xi3) < 1.0e-15;

    Vec3d x3 = v_min(x1, x2);
    success += x3[0] == (xi1 < xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 < xii2 ? xi2 : xii2)
               && x3[2] == (xi3 < xii3 ? xi3 : xii3);

    x3 = v_max(x1, x2);
    success += x3[0] == (xi1 > xii1 ? xi1 : xii1) 
               && x3[1] == (xi2 > xii2 ? xi2 : xii2)
               && x3[2] == (xi3 > xii3 ? xi3 : xii3);

    Vec3d x4;
    orthogonal(x2, x3, x4);
    success += abs(dot(x2, x3)) < 1.0e-15;

    x3 = cross(x1, x2);
    success += abs(x3[0] - xi2*xii3 + xi3*xii2) < 1.0e-15
               && abs(x3[1] - xi3*xii1 + xi1*xii3) < 1.0e-15
               && abs(x3[2] - xi1*xii2 + xi2*xii1) < 1.0e-15;

    double theta, phi, r;
    x1.get_spherical(theta, phi, r);
    success += abs(theta - acos(xi3/sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3))) < 1.0e-15
               && abs(phi - atan(xi2/xi1)) < 1.0e-15 
               && abs(r - sqrt(xi1*xi1 + xi2*xi2 + xi3*xi3)) < 1.0e-15;

    x3.set_spherical(theta, phi, r);
    success += abs(x3[0] - xi1) < 1.0e-15 
               && abs(x3[1] - xi2) < 1.0e-15 
               && abs(x3[2] - xi3) < 1.0e-15;
  }
  if(success != 8)
  {
    cout << "Failure in test of Vec3d Vector operations";
    return 1;
  }
  success = 0;

  ////////////////////////////////////////////////
  //              F I N A L I Z E
  ////////////////////////////////////////////////

  cout << "Performance time: " << t.get_secs();
  cout << endl << "ArithVec and derived vector classes have been tested successfully";

  return 0;
}
