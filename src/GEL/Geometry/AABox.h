/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file AABox.h
 * @brief Axis aligned bounding box class.
 */

#ifndef __GEOMETRY_AABOX__H
#define __GEOMETRY_AABOX__H

#include <iostream>
#include <vector>
#include "Triangle.h"

namespace Geometry
{
  const float DIST_THRESH = 5.0e-4f;

  class AABox
  {
	CGLA::Vec3f pmin, pmax, interior_point;
  public:

	AABox() {}

	AABox(const CGLA::Vec3f& _pmin, const CGLA::Vec3f& _pmax,
		  const CGLA::Vec3f& _interior_point):

	pmin(_pmin), pmax(_pmax), interior_point(_interior_point)
	{
      for(int i=0;i<3;++i)
		if((pmax[i]-pmin[i]) < DIST_THRESH)
		{
			pmax[i] += DIST_THRESH/2.0f;
			pmin[i] -= DIST_THRESH/2.0f;
		}
	  assert(pmin.all_le(interior_point));
	  assert(pmax.all_ge(interior_point));
	}

	const CGLA::Vec3f& get_pmin() const {return pmin;}

	const CGLA::Vec3f& get_pmax() const {return pmax;}

	bool intersect(const CGLA::Vec3f&, const CGLA::Vec3f&) const;

	void minmax_sq_dist(const CGLA::Vec3f& p, float& dmin, float& dmax) const;

	static AABox box_triangle(const Triangle&);

	static AABox box_and_split(const std::vector<Triangle>& invec,
							   std::vector<Triangle>& lvec,
							   std::vector<Triangle>& rvec);
														 
  };
}

#endif
