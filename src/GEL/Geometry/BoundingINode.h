/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BoundingINode.h
 * @brief Interior node class of a bounding hierarchy
 */

#ifndef __GEOMETRY_BOUNDINGINODE_H
#define __GEOMETRY_BOUNDINGINODE_H

#include "Ray.h"
#include "BoundingNode.h"

namespace Geometry
{

/** Interior node of bounding box tree. Contains
		pointers to left and right tree. */
template<class BoxType>
class BoundingINode: public BoundingNode<BoxType>
{
	BoundingNode<BoxType>* left;
	BoundingNode<BoxType>* right;
 public:

	BoundingINode(const BoxType& box,
								BoundingNode<BoxType>* _left,
								BoundingNode<BoxType>* _right): 
		BoundingNode<BoxType>(box), left(_left), right(_right) {}

	virtual ~BoundingINode() {delete left; delete right;}

	bool intersect(const CGLA::Vec3f&,const CGLA::Vec3f&,float&) const; 
	void intersect(Ray&  r) const; 
	int intersect_cnt(const CGLA::Vec3f&,const CGLA::Vec3f&) const;
 
	const BoundingNode<BoxType>* get_left() const {return left;}
	const BoundingNode<BoxType>* get_right() const {return right;}
};

}
#endif
