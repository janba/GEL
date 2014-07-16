/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "BoundingNode.h"
#include "BoundingINode.h"
#include "BoundingLNode.h"


using namespace std;
using namespace CGLA;

namespace Geometry
{

template<class BoxType>
void BoundingNode<BoxType>::sq_distance(const Vec3f& p, 
										float& dmin, float& dmax,
										float& s) const
{
  BoxType::minmax_sq_dist(p, dmin, dmax);
	s = 0;
}

template<class BoxType>
BoundingNode<BoxType>* 
BoundingNode<BoxType>::build(std::vector<Triangle>& triangles)
{
	int N = triangles.size();
	if(N==1)
		{
			const Triangle& t = triangles[0];
#if USE_LEAF_BOXES
			return new BoundingLNode<BoxType>(BoxType::box_triangle(t), t);
#else
			return new BoundingLNode<BoxType>(t);
#endif

		}
	else
		{
			std::vector<Triangle> triangles_left;
			std::vector<Triangle> triangles_right;

			BoxType box = 
				BoxType::box_and_split(triangles, triangles_left, triangles_right);

			BoundingNode* left  = build(triangles_left);
			BoundingNode* right = build(triangles_right);

			BoundingNode<BoxType>* bn = new BoundingINode<BoxType>(box, left, right);
			return bn;
		}
}

template class BoundingNode<AABox>;
/*
template BoundingNode<AABox>* 
BoundingNode<AABox>::build(std::vector<Triangle>& triangles);
*/
template class BoundingNode<OBox>;
/*template BoundingNode<OBox>* 
BoundingNode<OBox>::build(std::vector<Triangle>& triangles);
*/
}
