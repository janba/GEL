/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BoundingNode.h
 * @brief Abstract ancestor of the nodes of a bounding hierarchy.
 */

#ifndef __GEOMETRY_BOUNDINGNODE_H
#define __GEOMETRY_BOUNDINGNODE_H

#include <vector>
#include "../CGLA/Vec3f.h"
#include "Ray.h"
#include "Triangle.h"
#include "AABox.h"
#include "OBox.h"

namespace Geometry
{

/// Abstract BOUNDINGNODE node.
template<class BoxType>
class BoundingNode: public BoxType
{
 public:
	
	BoundingNode(const BoxType& box): BoxType(box) {}
	virtual ~BoundingNode() {}

	/// Count number of intersections from a point in a given direction
	virtual int intersect_cnt(const CGLA::Vec3f&,const CGLA::Vec3f&) const=0;

	/// Find the surface intersection point
	virtual bool intersect(const CGLA::Vec3f&,const CGLA::Vec3f&,float&) const=0;

	virtual void intersect(Ray& r) const = 0;

	/** For a given point, return the min and max square distance and the
			sign. Non-leafs return zero for the sign. Leaves
			return the square of the true distance as both min and max. */
	virtual void sq_distance(const CGLA::Vec3f&, float&, float&, float&) const;
	
	static BoundingNode* build(std::vector<Triangle>& triangles);
};

}
#endif
