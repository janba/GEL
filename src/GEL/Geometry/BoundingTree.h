/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BoundingTree.h
 * @brief Template representing a bounding hierarchy.
 */

#ifndef __GEOMETRY_BOUNDINGTREE_H
#define __GEOMETRY_BOUNDINGTREE_H

#include "BoundingNode.h"
#include "BoundingLNode.h"
#include "BoundingINode.h"
#include "Ray.h"

namespace Geometry
{

/** Template representing a bounding hierarchy. The argument should be the bounding box type - 
    either ABox or OOBox */
template<class BoxType>
class BoundingTree
{
 public:

	typedef BoundingNode<BoxType> Node;
	typedef BoundingLNode<BoxType> LeafNode;
	typedef BoundingINode<BoxType> IntNode;
	
	Node* root;

 public:

	BoundingTree(): root(0) {}

	~BoundingTree() {delete root;}

	void build(std::vector<Triangle>& triangles);

	bool intersect(const CGLA::Vec3f&,const CGLA::Vec3f&,float&) const;

	void intersect(Ray& r) const;

	int intersect_cnt(const CGLA::Vec3f&,const CGLA::Vec3f&) const;

	float compute_signed_distance(const CGLA::Vec3f& p,	float=FLT_MAX) const;
};

}

#endif
