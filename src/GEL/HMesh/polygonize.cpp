//
//  polygonize.cpp
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/Implicit.h>
#include <GEL/Geometry/Neighbours.h>
#include <GEL/HMesh/polygonize.h>
#include <GEL/HMesh/mesh_optimization.h>
#include <GEL/HMesh/smooth.h>
#include <GEL/HMesh/cleanup.h>
#include <GEL/HMesh/triangulate.h>

using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace HMesh;

namespace {
    const Vec3d hex_faces[6][4] = {
        {Vec3d(-0.5,-0.5,-0.5), Vec3d(-0.5,0.5,-0.5),   Vec3d(-0.5,0.5,0.5),    Vec3d(-0.5,-0.5,0.5)},
        {Vec3d(0.5, 0.5,-0.5),  Vec3d(0.5,-0.5,-0.5),   Vec3d(0.5,-0.5,0.5),    Vec3d(0.5,0.5,0.5)},
        {Vec3d( 0.5,-0.5, -0.5),Vec3d(-0.5,-0.5, -0.5), Vec3d(-0.5,-0.5, 0.5),  Vec3d(0.5,-0.5, 0.5)},
        {Vec3d(-0.5,0.5, -0.5), Vec3d(0.5,0.5, -0.5),   Vec3d(0.5,0.5, 0.5),    Vec3d(-0.5,0.5, 0.5)},
        {Vec3d(-0.5,-0.5,-0.5), Vec3d(0.5,-0.5,-0.5),   Vec3d(0.5,0.5,-0.5),    Vec3d(-0.5,0.5,-0.5)},
        {Vec3d( 0.5,-0.5,0.5),  Vec3d(-0.5,-0.5,0.5),   Vec3d(-0.5,0.5,0.5),    Vec3d(0.5,0.5,0.5)}
    };
}

namespace HMesh
{
    float clamp_interpolate(const RGrid<float>& grid, const CGLA::Vec3d& _v)
    {
        const Vec3i c0i = v_min(grid.get_dims()-Vec3i(2), v_max(Vec3i(_v), Vec3i(0)));
        const Vec3d v = v_min(Vec3d(grid.get_dims()-Vec3i(1)), v_max(_v,Vec3d(0)));
        
        const float alpha = v[0] - float(c0i[0]);
        const float beta  = v[1] - float(c0i[1]);
        const float gamm  = v[2] - float(c0i[2]);
        float m_alpha = 1.0 - alpha;
        float m_beta  = 1.0 - beta;
        float m_gamm  = 1.0 - gamm;
        array<float, 8> weights = {
            m_alpha*m_beta*m_gamm,
            alpha*m_beta*m_gamm,
            m_alpha*beta*m_gamm,
            alpha*beta*m_gamm,
            m_alpha*m_beta*gamm,
            alpha*m_beta*gamm,
            m_alpha*beta*gamm,
            alpha*beta*gamm};
        
        float f = 0;
        for(int i=0;i<8;++i)
            f += weights[i]*grid[c0i+Geometry::CubeCorners8i[i]];
        
        return f;
    }

    Vec3f clamp_trilin_grad(const RGrid<float>& grid, const CGLA::Vec3d& _v)
    {
        const Vec3i c0i = v_min(grid.get_dims()-Vec3i(2), v_max(Vec3i(_v), Vec3i(0)));
        const Vec3d v = v_min(Vec3d(grid.get_dims()-Vec3i(1)), v_max(_v,Vec3d(0)));

        const float alpha = v[0] - float(c0i[0]);
        const float beta  = v[1] - float(c0i[1]);
        const float gamm  = v[2] - float(c0i[2]);
        float m_alpha = 1.0 - alpha;
        float m_beta  = 1.0 - beta;
        float m_gamm  = 1.0 - gamm;
        array<float, 8> dxweights = {
            -m_beta*m_gamm,
            m_beta*m_gamm,
            -beta*m_gamm,
            beta*m_gamm,
            -m_beta*gamm,
            m_beta*gamm,
            -beta*gamm,
            beta*gamm};

        array<float, 8> dyweights = {
            -m_alpha*m_gamm,
            -alpha*m_gamm,
            m_alpha*m_gamm,
            alpha*m_gamm,
            -m_alpha*gamm,
            -alpha*gamm,
            m_alpha*gamm,
            alpha*gamm};

        array<float, 8> dzweights = {
            -m_alpha*m_beta,
            -alpha*m_beta,
            -m_alpha*beta,
            -alpha*beta,
            m_alpha*m_beta,
            alpha*m_beta,
            m_alpha*beta,
            alpha*beta};

        Vec3f gf(0);
        for(int i=0;i<8;++i) {
            gf[0] += dxweights[i]*grid[c0i+Geometry::CubeCorners8i[i]];
            gf[1] += dyweights[i]*grid[c0i+Geometry::CubeCorners8i[i]];
            gf[2] += dzweights[i]*grid[c0i+Geometry::CubeCorners8i[i]];
        }
        return gf;
    }

      
    void polygonize(const RGrid<float>& grid, std::vector<CGLA::Vec3d>& quad_vertices,
                    float tau, bool high_is_inside)
    {
        auto is_inside = [&](const Vec3i& pi) {
            float val = grid[pi];
            return !std::isnan(val) && (high_is_inside == (val > tau));
        };
        auto is_outside = [&](const Vec3i& pi) {
            if (grid.in_domain(pi)) {
                float val = grid[pi];
                return !std::isnan(val) && (high_is_inside == (val <= tau));
            }
            return true;

        };
        
        quad_vertices.clear();
        for(Vec3i pi: Range3D(grid.get_dims())) {
            if(is_inside(pi)) {
                Vec3d p(pi);
                for (int nbr_idx = 0; nbr_idx < 6 ; ++ nbr_idx) {
                    Vec3i pni = pi + N6i[nbr_idx];
                    if(is_outside(pni))
                        for(int n=0;n<4;++n)
                            quad_vertices.push_back(p + hex_faces[nbr_idx][3-n]);
                }
            }
        }
    }

    void volume_polygonize(const XForm& xform, const Geometry::RGrid<float>& grid,
                           HMesh::Manifold& mani, float tau, bool make_triangles, bool high_is_inside)
    {
        const double delta = sqrt(3.0)/2.0;

        mani.clear();
        vector<Vec3d> quad_vertices;
        polygonize(grid, quad_vertices, tau, high_is_inside);
        vector<int> indices;
        vector<int> faces(quad_vertices.size()/4,4);
        for(int i=0;i<quad_vertices.size();++i)
            indices.push_back(i);
        build(mani, quad_vertices.size(),
                   quad_vertices[0].get(),
                   faces.size(),
                   &faces[0],
                   &indices[0]);
        
        stitch_mesh(mani, 0.001 * delta);

        if(make_triangles)
            triangulate(mani);
        mani.cleanup();
        
        struct Corner {
            float val;
            bool in_domain;
            Vec3d pos;
        };
        array<Corner,8> corners;
        for(auto v: mani.vertices()) {
            const Vec3i pi(mani.pos(v));
            int cnt = 0;
            float va=0;
            Vec3d pa(0.0);
            for (int i=0;i<8;++i) {
                Vec3i pbi = pi + CubeCorners8i[i];
                if (grid.in_domain(pbi)) {
                    corners[i].in_domain = true;
                    va += corners[i].val = grid[pbi];
                    pa += corners[i].pos = Vec3d(pbi);
                    ++cnt;
                }
                else corners[i].in_domain = false;
            }
            pa /= cnt;
            va /= cnt;
            
            Vec3d p(0.0);
            Vec3d n = normal(mani, v);
            cnt = 0;
            for (int i=0;i<8;++i)
                if (corners[i].in_domain) {
                    float vb = corners[i].val;
                    if(max(va,vb)>tau && min(va,vb)<=tau) {
                        float delta = vb - va;
                        Vec3d pb = corners[i].pos;
                        if(dot(n, pb-pa)>0.0) {
                            p += pa * (vb-tau)/delta - pb * (va-tau)/delta;
                            ++cnt;
                        }
                    }
                }
            mani.pos(v) = (cnt>0)? p / cnt : pa;

        }
        for(auto v: mani.vertices())
            mani.pos(v) = xform.inverse(mani.pos(v));
    }
    
}
