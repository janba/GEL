/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/Geometry/AABox.h>
#include <GEL/Geometry/BoundingINode.h>
#include <GEL/Geometry/OBox.h>

using namespace std;
using namespace CGLA;

namespace Geometry
{

template<class BoxType>
int BoundingINode<BoxType>::intersect_cnt(const CGLA::Vec3f& p , 
																 const CGLA::Vec3f& dir) const
{
	if(!BoxType::intersect(p,dir))
		return 0;

	int li = left->intersect_cnt(p,dir);
	int ri = right->intersect_cnt(p,dir);
	return li + ri;
}

template<class BoxType>
void BoundingINode<BoxType>::intersect(Ray& r) const 
{
	if(BoxType::intersect(r.origin,r.direction))
	{
			left->intersect(r);
			right->intersect(r);
	}
}

template<class BoxType>
bool BoundingINode<BoxType>::intersect(const CGLA::Vec3f& p , const CGLA::Vec3f& dir,
															float& tmin) const 
{
	if(!BoxType::intersect(p,dir))
		return false;

	float tminl= FLT_MAX, tminr= FLT_MAX;
	bool li = left->intersect(p,dir,tminl) && tminl > 0;
	bool ri = right->intersect(p,dir,tminr) && tminr > 0;
	if(li ||ri)
		{
			tminl = tminl>0 ? tminl : FLT_MAX;
			tminr = tminr>0 ? tminr : FLT_MAX;
			tmin = min(tminl,tminr);
			return true;
		}
	return false;
}

template class BoundingINode<AABox>;
template class BoundingINode<OBox>;

}
