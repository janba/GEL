/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "statistics.h"

#include "Mat2x2f.h"
#include "Mat3x3f.h"
#include "Mat4x4f.h"
#include "Mat2x2d.h"
#include "Mat3x3d.h"
#include "Mat4x4d.h"

using namespace std;

namespace CGLA
{
	template<class VT, class MT>
	VT covariance(const vector<VT>& vec, MT& C_out)
	{
		VT m = mean(vec);
		size_t n = vec.size();
		
		MT C(0);
		for(size_t i=0;i<n;++i)
		{
			MT B;
			VT v = vec[i]-m;
			outer_product(v,v,B);
			
			C += B;
		}
		C_out = C;
		
		return m;
	}
	
	template
	Vec2f covariance<Vec2f,Mat2x2f>(const vector<Vec2f>& vec, Mat2x2f& C_out);
	
	template 
	Vec3f covariance<Vec3f,Mat3x3f>(const vector<Vec3f>& vec, Mat3x3f& C_out);
	
	template 
	Vec4f covariance<Vec4f,Mat4x4f>(const vector<Vec4f>& vec, Mat4x4f& C_out);
	
	template 
	Vec2d covariance<Vec2d,Mat2x2d>(const vector<Vec2d>& vec, Mat2x2d& C_out);
	
	template 
	Vec3d covariance<Vec3d,Mat3x3d>(const vector<Vec3d>& vec, Mat3x3d& C_out);
	
	template 
	Vec4d covariance<Vec4d,Mat4x4d>(const vector<Vec4d>& vec, Mat4x4d& C_out);
}
