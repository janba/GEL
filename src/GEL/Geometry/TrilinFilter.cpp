/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "GradientFilter.h"
#include "TrilinFilter.h"
#include "Neighbours.h"


using namespace CGLA;

namespace Geometry
{
	template<class GridT>
	bool TrilinFilter<GridT>::in_domain(const CGLA::Vec3f& v) const
	{
		Vec3i c0i(v);
		return grid->in_domain(c0i) && grid->in_domain(c0i+Vec3i(1));
	}
    
	template<class GridT>
	float TrilinFilter<GridT>::operator()(const CGLA::Vec3f& v) const
	{
		Vec3i c0i(v);
        
		const float alpha = v[0] - float(c0i[0]);
		const float beta  = v[1] - float(c0i[1]);
		const float gamm  = v[2] - float(c0i[2]);
        
		float m_alpha = 1.0 - alpha;
		float m_beta  = 1.0 - beta;
		float m_gamm  = 1.0 - gamm;
        
		float voxels[8];
		for(int i=0;i<8;++i)
			voxels[i] = (*grid)[c0i+Geometry::CubeCorners8i[i]];
        
		float f = 0.0f;
		f += (m_alpha*m_beta*m_gamm)*voxels[0];
		f += (alpha*m_beta*m_gamm)*voxels[1];
		f += (m_alpha*beta*m_gamm)*voxels[2];
		f += (alpha*beta*m_gamm)*voxels[3];
		f += (m_alpha*m_beta*gamm)*voxels[4];
		f += (alpha*m_beta*gamm)*voxels[5];
		f += (m_alpha*beta*gamm)*voxels[6];
		f += (alpha*beta*gamm)*voxels[7];
		
		return f;
	}
    
	template<class GridT>
	Vec3f TrilinFilter<GridT>::grad(const CGLA::Vec3f& v) const
	{
		Vec3i c0i(v);
        
		const float alpha = v[0] - float(c0i[0]);
		const float beta  = v[1] - float(c0i[1]);
		const float gamm  = v[2] - float(c0i[2]);
        
		float m_alpha = 1.0 - alpha;
		float m_beta  = 1.0 - beta;
		float m_gamm  = 1.0 - gamm;
        
		float voxels[8];
		for(int i=0;i<8;++i)
			voxels[i] = (*grid)[c0i+Geometry::CubeCorners8i[i]];
        
		float fx = 0.0f;
		fx += -m_beta*m_gamm*float(voxels[0]);
		fx += m_beta*m_gamm*float(voxels[1]);
		fx += -beta*m_gamm*float(voxels[2]);
		fx += beta*m_gamm*float(voxels[3]);
		fx += -m_beta*gamm*float(voxels[4]);
		fx += m_beta*gamm*float(voxels[5]);
		fx += -beta*gamm*float(voxels[6]);
		fx += beta*gamm*float(voxels[7]);
        
		float fy = 0.0f;
		fy += -m_alpha*m_gamm*float(voxels[0]);
		fy += -alpha*m_gamm*float(voxels[1]);
		fy += m_alpha*m_gamm*float(voxels[2]);
		fy += alpha*m_gamm*float(voxels[3]);
		fy += -m_alpha*gamm*float(voxels[4]);
		fy += -alpha*gamm*float(voxels[5]);
		fy += m_alpha*gamm*float(voxels[6]);
		fy += alpha*gamm*float(voxels[7]);
        
		float fz = 0.0f;
		fz += -m_alpha*m_beta*float(voxels[0]);
		fz += -alpha*m_beta*float(voxels[1]);
		fz += -m_alpha*beta*float(voxels[2]);
		fz += -alpha*beta*float(voxels[3]);
		fz += m_alpha*m_beta*float(voxels[4]);
		fz += alpha*m_beta*float(voxels[5]);
		fz += m_alpha*beta*float(voxels[6]);
		fz += alpha*beta*float(voxels[7]);
        
		return Vec3f(fx, fy, fz);
	}
    
    template class  TrilinFilter<RGrid<float> >;
    template class  TrilinFilter<RGrid<unsigned char> >;
    template class  TrilinFilter<RGrid<unsigned short> >;
    template class  TrilinFilter<RGrid<short> >;
}
