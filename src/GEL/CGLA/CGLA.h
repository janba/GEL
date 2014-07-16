/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file CGLA.h
 * @brief CGLA main header: contains a number of constants and function decs.
 */

#ifndef __CGLA_CGLA_H__
#define __CGLA_CGLA_H__

#if (_MSC_VER >= 1200)
#pragma warning (disable: 4244 4800)
#endif

#include <cmath>
#include <cfloat>
#include <climits>
#include <cassert>
#include <algorithm>
#include <functional>

#ifdef _MSC_VER // if visual studio
#include <float.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#endif

#ifndef DEGREES_TO_RADIANS
#define DEGREES_TO_RADIANS (M_PI / 180.0)
#endif

namespace CGLA 
{
		inline float cgla_nan() 
		{
				static const float cgla_nan_value = log(-1.0f);
				return cgla_nan_value;
		}
		
  /** Procedural definition of NAN */  
#define CGLA_NAN cgla_nan()

  /** NAN is used for initialization of vectors and matrices
      in debug mode */
#define CGLA_INIT_VALUE cgla_nan()

  /** Numerical constant representing something large.
      value is a bit arbitrary */
  const double BIG=10e+30;

  /** Numerical constant represents something extremely small.
      value is a bit arbitrary */
  const double MINUTE=10e-30;

  /** Numerical constant represents something very small.
      value is a bit arbitrary */
  const double TINY=3e-7;
	
  /** Numerical constant represents something small.
      value is a bit arbitrary */
  const double SMALL=10e-2;

  /** The GEL pseudo-random number generator uses
      UINT_MAX as RAND_MAX to avoid mod operations. */
  const unsigned int GEL_RAND_MAX=UINT_MAX;

		inline double sqrt3()
		{
				static const double sqrt3_val = sqrt(3.0);
				return sqrt3_val;
		}

#define SQRT3 sqrt3()

  /// Useful enum that represents coordiante axes.
  enum Axis {XAXIS=0,YAXIS=1,ZAXIS=2};

  inline bool isnan(double x) {
#ifdef _MSC_VER // if visual studio
      return _isnan(x);
#else
      return std::isnan(x);
#endif
  }

  /// Template for a function that squares the argument.
  template <class Scalar>
  inline Scalar sqr(Scalar x) {///
    return x*x;}
	
  /// Scalar template for a function that returns the cube of the argument.
  template <class Scalar>
  inline Scalar qbe(Scalar x) {///
    return x*x*x;}

  template <class Scalar>
  inline bool is_zero(Scalar x)	{return (x > -MINUTE && x < MINUTE);}

  template <class Scalar>
  inline bool is_tiny(Scalar x)	{return (x > -TINY && x < TINY);}

  /** What power of 2 ?. if x is the argument, find the largest 
      y so that 2^y <= x */
  inline int two_to_what_power(unsigned int x) 
  {
    if (x<1) 
      return -1;
    int i = 0;
    while (x != 1) {x>>=1;i++;}
    return i;
  }

#ifdef __sgi
  inline int round(float x) {return int(rint(x));}
#else
  inline int round(float x) {return int(x+0.5);}
#endif

  template<class T>
  inline T sign(T x) {return x>=T(0) ? 1 : -1;}

  /// Integer power function with O(log(n)) complexity
  template<class T>
  inline T int_pow(T a, unsigned int n)
  {
    T result = static_cast<T>(1);
    for(; n > 0; n >>= 1)
    {
      if(n & 1) result = result*a;
      a *= a;
    }
    return result;
  }
  
  /// Function that seeds the GEL pseudo-random number generator
  void gel_srand(unsigned int seed);

  /** GEL provides a linear congruential pseudo-random number 
      generator which is optimized for speed. This version allows 
      an integer argument which is useful for grid-based noise
      functions. */
  unsigned int gel_rand(unsigned int k);

  /** GEL provides a linear congruential pseudo-random number 
      generator which is optimized for speed. This means
      that GEL_RAND_MAX==UINT_MAX. */
  unsigned int gel_rand();

	/** raw_assign takes a CGLA vector, matrix or whatever has a get() function
			as its first argument and a raw pointer to a (presumed scalar) entity 
			as the second argument. the contents dereferenced by the pointer is 
			copied to the entity given as first argument. */
  template<class T, class S>
  void raw_assign(T& a,  const S* b)
  {
    memcpy(static_cast<void*>(a.get()),static_cast<const void*>(b),sizeof(T));
  }
	
}

#endif
