/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BoundingLNode.h
 * @brief Leaf node of a bounding hierarchy.
 */

#ifndef __GEOMETRY_BOUNDINGLNODE_H
#define __GEOMETRY_BOUNDINGLNODE_H

#include "Ray.h"
#include "BoundingNode.h"

#define USE_LEAF_BOXES 0

namespace Geometry
{

/// Leaf node of a bounding box tree. References triangle.
template<class BoxType>
class BoundingLNode: public BoundingNode<BoxType>
{
	Triangle tri;
 public:

#if USE_LEAF_BOXES
	BoundingLNode(const BoxType& box,
								const Triangle& _tri): BoundingNode<BoxType>(box), tri(_tri) {}
#else
	BoundingLNode(const Triangle& _tri): BoundingNode<BoxType>(BoxType()), tri(_tri) {}
#endif

	virtual ~BoundingLNode() {}

	bool intersect(const CGLA::Vec3f&,const CGLA::Vec3f&,float&) const;
	void intersect(Ray&r) const;
	int intersect_cnt(const CGLA::Vec3f&,const CGLA::Vec3f&) const;

	void sq_distance(const CGLA::Vec3f&, float&, float&, float&) const;

	virtual bool is_leaf() const {return true;}

	const Triangle& get_tri() const {return tri;}
};

template<class BoxType>
inline bool BoundingLNode<BoxType>::intersect(const CGLA::Vec3f& p, 
																							const CGLA::Vec3f& dir,
																							float& tmin) const
{
#if USE_LEAF_BOXES
	if(BoxType::intersect(p,dir))
 		return tri.intersect(p,dir,tmin);
 	return false;
#else
	return tri.intersect(p,dir,tmin);
#endif

}

template<class BoxType>
inline void BoundingLNode<BoxType>::intersect(Ray& r) const
{
		CGLA::Vec3f p = r.origin;
		CGLA::Vec3f d = r.direction;
		float t = FLT_MAX;
		if(tri.intersect(p,d,t) && t < r.dist)
		{
				r.has_hit = true;
				r.dist = t;
				r.hit_pos = p + t*d;
				r.hit_normal = tri.get_face_norm();
		}
}

template<class BoxType>
inline int BoundingLNode<BoxType>::intersect_cnt(const CGLA::Vec3f& p, 
																								 const CGLA::Vec3f& dir) const
{
#if USE_LEAF_BOXES
	float tmin=1.0e33f;
	if(BoxType::intersect(p,dir) && 
		 tri.intersect(p,dir,tmin) &&
		 tmin > 0)
		return 1;
	return 0;
#else
	float tmin=1e33f;
	if(tri.intersect(p,dir,tmin) &&
		 tmin > 0)
		return 1;
	return 0;
#endif
}

template<class BoxType>
inline void BoundingLNode<BoxType>::sq_distance(const CGLA::Vec3f& p, 
																								float& dmin, 
																								float& dmax, float& s) const
{
 	bool did_work = tri.signed_distance(p,dmax,s);
	if(!did_work) std::cout << dmax << std::endl;
	assert(did_work);
 	dmin = dmax; 
}

}
#endif
