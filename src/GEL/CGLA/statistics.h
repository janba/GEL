/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file statistics.h
 * Compute mean and covariance of CGLA vectors - simple multivariate statistics.
 */

#ifndef __CGLA_STATISTICS_H__
#define __CGLA_STATISTICS_H__

#if (_MSC_VER >= 1200)
#pragma warning (disable: 4018 4244 4800)
#endif

#include <vector>

namespace CGLA
{
		template<class VT>
				VT mean(const std::vector<VT>& vec)
				{
						VT v(0);
						for(unsigned int i=0;i<vec.size();++i)
								v += vec[i];
						v /= vec.size();

						return v;
				}


		/** Function that computes the covariance of a set of points.
				This function returns the mean, and, upon completion, the
				final argument contains the covariance matrix.

				This template is instantiated for Vec3f, Vec2f, and Vec4f. */
			
		template<class VT, class MT>
				VT covariance(const std::vector<VT>& vec, MT& C_out);
}




#endif
