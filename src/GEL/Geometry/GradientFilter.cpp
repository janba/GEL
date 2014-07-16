/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */


//#include <iostream>
#include "GradientFilter.h"
#include "../Geometry/Neighbours.h"
#include "../Geometry/RGrid.h"

using namespace std;
using namespace CGLA;
using namespace Geometry;


namespace Geometry
{
	template<class GridT>
	bool GradientFilter<GridT>::in_domain(const CGLA::Vec3i& p) const
	{
		const Vec3i& dims = grid->get_dims();
		for(int i=0; i<3; ++i)
			if ((p[i]-1)<0 || (p[i]+1) >= dims[i])
				return false;
		return true;
	}

	template<class GridT>
	bool GradientFilter<GridT>::map(const Vec3i& p, Vec3f& g) const
	{
//		bool pit, mit;
		
		float vxm = (*grid)[p+N6i[0]];
		float vxp = (*grid)[p+N6i[1]];
		float vym = (*grid)[p+N6i[2]];
		float vyp = (*grid)[p+N6i[3]];
		float vzm = (*grid)[p+N6i[4]];
		float vzp = (*grid)[p+N6i[5]];
		float dx = (float(vxp)-vxm)/2.0f;
		float dy = (float(vyp)-vym)/2.0f;
		float dz = (float(vzp)-vzm)/2.0f;
		g = Vec3f(dx,dy,dz);
		return true;
	}


	template class GradientFilter<Geometry::RGrid<float> >;
	template class GradientFilter<Geometry::RGrid<unsigned char> >;
	template class GradientFilter<Geometry::RGrid<unsigned short> >;
}
