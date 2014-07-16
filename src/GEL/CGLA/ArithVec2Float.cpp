/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>
#include <algorithm>
#include "Vec2f.h"
#include "Vec2d.h"
#include "Mat2x2f.h"
#include "Mat2x2d.h"
#include "CGLA.h"

namespace CGLA {
  using namespace std;

  template<class T, class V>
  bool linear_combine(const ArithVec2Float<T,V>& a,
		      const ArithVec2Float<T,V>& b,
		      const ArithVec2Float<T,V>& c,
		      T& x,
		      T& y) 
  {
    Vec2d xy = invert(Mat2x2d(a[0],b[0],a[1],b[1])) * Vec2d(c[0], c[1]);
    x = xy[0];
    y = xy[1];
    return true;
  }
	

  template class ArithVec2Float<float, Vec2f>;
  template class ArithVec2Float<double, Vec2d>;

  template bool 
  linear_combine<double,Vec2d>(const ArithVec2Float<double,Vec2d>&,
			       const ArithVec2Float<double,Vec2d>&,
			       const ArithVec2Float<double,Vec2d>&,
			       double&,
			       double&);

  template bool 
  linear_combine<float,Vec2f>(const ArithVec2Float<float,Vec2f>&,
			      const ArithVec2Float<float,Vec2f>&,
			      const ArithVec2Float<float,Vec2f>&,
			      float&,
			      float&);


}
