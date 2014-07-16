/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file OBox.h
 * @brief An oriented bounding box.
 */

#ifndef __GEOMETRY_OBOX__H
#define __GEOMETRY_OBOX__H

#include <iostream>
#include <vector>
#include "Triangle.h"
#include "AABox.h"

namespace Geometry
{

class OBox
{
	const CGLA::Mat3x3f R;
	const AABox aabox;

public:
	OBox() {}

	OBox(const CGLA::Mat3x3f& _R, const AABox& _aabox):
		R(_R), aabox(_aabox) {}

	bool intersect(const CGLA::Vec3f&, const CGLA::Vec3f&) const;

	void minmax_sq_dist(const CGLA::Vec3f& p, float& dmin, float& dmax) const
	{
		aabox.minmax_sq_dist(R * p, dmin, dmax);
	}
	

	static OBox box_triangle(const Triangle&);

	static OBox box_and_split(const std::vector<Triangle>& invec,
														 std::vector<Triangle>& lvec,
														 std::vector<Triangle>& rvec);
														 
	const CGLA::Mat3x3f& get_rotation() const { return R; }
	const AABox& get_aabox() const { return aabox; }

/* 	const CGLA::Vec3f& get_pmin() const {assert(0); return CGLA::Vec3f(0);} */

/* 	const CGLA::Vec3f& get_pmax() const {assert(0); return CGLA::Vec3f(0);} */

};

}
#endif
