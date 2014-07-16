/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <cfloat>
#include "AABox.h"

using namespace std;
using namespace CGLA;

namespace Geometry
{

bool AABox::intersect(const CGLA::Vec3f& p, const CGLA::Vec3f& dir) const
{
	Vec3f t0,t1;
	for(int i=0;i<3;++i)
		{
			t0[i] = (pmin[i]-p[i])/dir[i];
			t1[i] = (pmax[i]-p[i])/dir[i];
		}
	Vec3f tin = v_min(t0, t1);
	Vec3f tout = v_max(t0,t1);
	float tmin = max(tin[0], max(tin[1], tin[2]));
	float tmax = min(tout[0], min(tout[1], tout[2]));

	return ( (tmin-CGLA::TINY) < (tmax+CGLA::TINY));
}

void AABox::minmax_sq_dist(const CGLA::Vec3f& p, 
													 float& dmin, float& dmax) const
{
	const Vec3f a = 0.5*pmax-0.5*pmin;
	const Vec3f p0 = pmin + a;
 	Vec3f d = p-p0;
 	Vec3f f(d);

	for(int i=0;i<3;++i)
		{
 			if(f[i]>=0) 
 				f[i] = p[i]-pmin[i];
 			else
 				f[i] = p[i]-pmax[i];
			
			if(d[i]<-a[i])
				d[i] = p[i]-pmin[i];
			else if(d[i]>a[i])
				d[i] = p[i]-pmax[i];
			else 
				d[i] = 0;
		}
	dmin = sqr_length(d);
//  	dmax = sqr_length(p-interior_point);
	dmax = sqr_length(f);
	assert(dmin<=dmax);
}

AABox AABox::box_triangle(const Triangle& t)
{
	return AABox(t.get_pmin(), t.get_pmax(), t.get_centre());
}


AABox AABox::box_and_split(const std::vector<Triangle>& invec,
													 std::vector<Triangle>& lvec,
													 std::vector<Triangle>& rvec)
{
	const size_t N = invec.size();
	Vec3f tri_pmin(FLT_MAX), tri_pmax(-FLT_MAX);
			
	for(size_t i=0;i<N;++i)
		{
			tri_pmin = v_min(invec[i].get_pmin(), tri_pmin);
			tri_pmax = v_max(invec[i].get_pmax(), tri_pmax);
		}
	Vec3f diff = tri_pmax - tri_pmin;

	// Find the point closest to the centre.
	Vec3f centre = tri_pmin + diff;
	Vec3f centre_close = invec[0].get_v0();
	float min_dist = FLT_MAX;
	for(size_t i=0;i<N;++i)
		{
			Vec3f v0 = invec[i].get_v0();
			Vec3f v1 = invec[i].get_v1();
			Vec3f v2 = invec[i].get_v2();
			float sl0 = sqr_length(centre-v0);
			if(sl0 < min_dist)
				{
					min_dist = sl0;
					centre_close = v0;
				}
			float sl1 = sqr_length(centre-v1);
			if(sl1 < min_dist)
				{
					min_dist = sl1;
					centre_close = v1;
				}
			float sl2 = sqr_length(centre-v2);
			if(sl2 < min_dist)
				{
					min_dist = sl2;
					centre_close = v2;
				}
		}

	int k;
	if(diff[0]>diff[1])
		{
			if(diff[0]>diff[2]) 
				k = 0;
			else 
				k = 2;
		}
	else
		{
			if(diff[1]>diff[2]) 
				k = 1;
			else 
				k = 2;
		}

	float thresh = diff[k]/2.0f + tri_pmin[k];

 	for(size_t i=0;i<N;++i)
		{
			if(invec[i].get_centre()[k] > thresh)
				rvec.push_back(invec[i]);
			else
				lvec.push_back(invec[i]);
		}
	if(lvec.empty() || rvec.empty())
		{
			lvec.clear();
			lvec.insert(lvec.end(),
									invec.begin(),
									invec.begin()+N/2);
			rvec.clear();
			rvec.insert(rvec.end(),
									invec.begin()+N/2,
									invec.end());
		}
	assert(!lvec.empty());
	assert(!rvec.empty());
	assert(lvec.size()+rvec.size() == invec.size());
	return AABox(tri_pmin, tri_pmax, centre_close);
}

}
