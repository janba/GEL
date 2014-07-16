//
//  Implicit.cpp
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#include "../CGLA/CGLA.h"
#include "../Geometry/Neighbours.h"
#include "../Geometry/GridAlgorithm.h"
#include "XForm.h"
#include "Implicit.h"

using namespace std;
using namespace CGLA;

namespace Geometry
{
    void Implicit::push_to_surface(Vec3d& p, double tau, double max_dist) const
    {
        Vec3d g = grad(p);
        
        double sl = sqr_length(g);
        double d = (eval(p)-tau)/sl;
        Vec3d disp = g*d;
        double disp_len = length(disp);
        double clamped_disp_len = min(max_dist, disp_len);
        p = p - clamped_disp_len*disp/disp_len;
    }
    
    XForm grid_sample(const Implicit& imp, const CGLA::Vec3d& llf, const CGLA::Vec3d& urt,
                      Geometry::RGrid<float>& grid)
    {
        XForm xform(llf,urt, grid.get_dims(), 0.0);
        for(int i=0;i<xform.get_dims()[0];++i)
            for(int j=0;j<xform.get_dims()[1];++j)
                for(int k=0;k<xform.get_dims()[2];++k)
                {
                    Vec3d p = xform.inverse(Vec3d(i,j,k));
                    float f = imp.eval(p);
                    grid[Vec3i(i,j,k)] = f;
                }
        return xform;
    }
    
    
    float interpolate(const RGrid<float>& grid, const CGLA::Vec3d& _v)
    {
        Vec3d v = v_min(Vec3d(grid.get_dims()-Vec3i(1)), v_max(_v,Vec3d(0)));
        
        Vec3i c0i(v);
        
        const float alpha = v[0] - float(c0i[0]);
        const float beta  = v[1] - float(c0i[1]);
        const float gamm  = v[2] - float(c0i[2]);
        float m_alpha = 1.0 - alpha;
        float m_beta  = 1.0 - beta;
        float m_gamm  = 1.0 - gamm;
        float weights[8];
        weights[0] = (m_alpha*m_beta*m_gamm);
        weights[1] = (alpha*m_beta*m_gamm);
        weights[2] = (m_alpha*beta*m_gamm);
        weights[3] = (alpha*beta*m_gamm);
        weights[4] = (m_alpha*m_beta*gamm);
        weights[5] = (alpha*m_beta*gamm);
        weights[6] = (m_alpha*beta*gamm);
        weights[7] = (alpha*beta*gamm);
        
        float f = 0;
        for(int i=0;i<8;++i)
            if(weights[i]>1e-10)
                f += weights[i]*grid[c0i+Geometry::CubeCorners8i[i]];
        
        return f;
    }
    
    
    double VolumetricImplicit::eval(const CGLA::Vec3d& p) const
    {
        return interpolate(grid,xform.apply(p));
    }
    
    CGLA::Vec3d VolumetricImplicit::grad(const CGLA::Vec3d& _p) const
    {
        Vec3d p = xform.apply(_p);
        Vec3d g(interpolate(grid,p+Vec3d(1,0,0))-interpolate(grid,p-Vec3d(1,0,0)),
                interpolate(grid,p+Vec3d(0,1,0))-interpolate(grid,p-Vec3d(0,1,0)),
                interpolate(grid,p+Vec3d(0,0,1))-interpolate(grid,p-Vec3d(0,0,1)));
        g *= 0.5 / xform.inv_scale();
        return g;
    }
    
}