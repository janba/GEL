/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file QEM.h
 * @brief Quadric Error Metrics. To be used when simplifying.
 */

#ifndef __GEOMETRY_QEM_H
#define __GEOMETRY_QEM_H

#include <cfloat>
#include "../CGLA/Vec3d.h"
#include "../CGLA/Mat3x3d.h"


namespace
{
	inline const CGLA::Mat3x3d direct_product(const CGLA::Vec3d& v0, const CGLA::Vec3d& v1)
		{
			CGLA::Mat3x3d m;
			for(int i=0;i<3;++i)
				for(int j=0;j<3;++j)
					m[i][j] = v0[i]*v1[j];
			return m;
		}
}

namespace Geometry
{
	class QEM
		{
			CGLA::Mat3x3d A;
			CGLA::Vec3d   b;
			double   c;
		public:
			
			QEM(): A(0), b(0), c(0) {}
			
			QEM(const CGLA::Vec3d& p0, const CGLA::Vec3d& n0, double w=1.0f):
				A(direct_product(n0,n0) * w), 
				b(-2*n0*dot(n0,p0) * w), 
				c(dot(p0,n0)*dot(p0,n0) * w) {}

			
			void operator+=(const QEM& q)
				{
					A += q.A;
					b += q.b;
					c += q.c;
				}
			
			float error(const CGLA::Vec3d& p) const
				{
					return dot(p,A*p) + dot(b,p)+ c;
				}

			double determinant() const
				{
					return CGLA::determinant(A);
				}
			
			const CGLA::Vec3d grad(const CGLA::Vec3d& p) const
				{
					return CGLA::Vec3d(2*A*p+b);
				}
			
			CGLA::Vec3d opt_pos(double QEM_thresh = 0.005, const CGLA::Vec3d& p0 = CGLA::Vec3d(0.0)) const;
			
		};
}

namespace GEO = Geometry;

#endif
