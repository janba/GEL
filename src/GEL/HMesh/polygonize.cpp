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
Vec3d hex_faces[6][4] = {{Vec3d(0,-0.5,-0.5),Vec3d(0,0.5,-0.5),Vec3d(0,0.5,0.5),Vec3d(0,-0.5,0.5)},
    {Vec3d(0, 0.5,-0.5),Vec3d(0,-0.5,-0.5),Vec3d(0,-0.5,0.5),Vec3d(0,0.5,0.5)},
    {Vec3d( 0.5,0, -0.5),Vec3d(-0.5,0, -0.5),Vec3d(-0.5,0, 0.5),Vec3d(0.5,0, 0.5)},
    {Vec3d(-0.5,0, -0.5),Vec3d(0.5,0, -0.5),Vec3d(0.5,0, 0.5),Vec3d(-0.5,0, 0.5)},
    {Vec3d(-0.5,-0.5,0),Vec3d(0.5,-0.5,0),Vec3d(0.5,0.5,0),Vec3d(-0.5,0.5,0)},
    {Vec3d( 0.5,-0.5,0),Vec3d(-0.5,-0.5,0),Vec3d(-0.5,0.5,0),Vec3d(0.5,0.5,0)}};
}

namespace HMesh
{
    float clamp_interpolate(const RGrid<float>& grid, CGLA::Vec3d& v)
    {
        v = v_min(Vec3d(grid.get_dims()-Vec3i(1)), v_max(v,Vec3d(0)));
        
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

    
    void polygonize(const RGrid<float>& grid, std::vector<CGLA::Vec3d>& quad_vertices,
                    float tau, bool high_is_inside)
    {
        auto is_inside = [&](const Vec3i& pi) {
            float val = grid[pi];
            return !isnan(val) && (high_is_inside == (val > tau));
        };
        auto is_outside = [&](const Vec3i& pi) {
            if (grid.in_domain(pi)) {
                float val = grid[pi];
                return isnan(val) || (high_is_inside == (val <= tau));
            }
            return true;

        };
        
        quad_vertices.clear();
        for(Vec3i pi: Range3D(grid.get_dims())) {
            if(is_inside(pi)) {
                Vec3d p(pi);
                for (int nbr_idx = 0; nbr_idx < 6 ; ++ nbr_idx) {
                    Vec3i pni = pi + N6i[nbr_idx];
                    if(is_outside(pni)) {
                        Vec3d pn = p + 0.5 * N6d[nbr_idx];
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(pn + hex_faces[nbr_idx][3-n]);
                    }
                }
            }
        }
    }


    bool towards_intersection(const Geometry::RGrid<float>& grid, Vec3d& p, const Vec3d& _p1, double tau) {
        Vec3d p1 = _p1;
        auto v0 = clamp_interpolate(grid, p);
        auto v1 = clamp_interpolate(grid, p1);
        double D = v1 - v0;
        
        if(abs(D)>1e-30) { // If the value at the two points are the same, we don't move
            double d = tau - v0;
            double r = d/D;
            if(r>1.0) // If the intersection is further out than p1, don't move.
                return false;
            if(r>=0.0) { // If the intersection is betweeen p0 and p1, move
                p = r * p1 + (1.0-r) * p;
                return true;
            }
            // Otherwise, there is no intersection.
        }
        return false;
    }
    
    
    void volume_polygonize(const XForm& xform, const Geometry::RGrid<float>& grid,
                           HMesh::Manifold& mani, float tau, bool make_triangles, bool high_is_inside)
    {
        const double delta = sqrt(3.0);

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

        for(int iter=0;iter<5; ++iter) {
            laplacian_smooth(mani,0.5);
            for(auto v: mani.vertices()) {
                Vec3d& p = mani.pos(v);
                const Vec3d n = normal(mani,v);
                const Vec3d p1 = p + delta * n;
                if(!towards_intersection(grid, p, p1, tau)) {
                    const Vec3d p1 = p - delta * n;
                    towards_intersection(grid, p, p1, tau);
                }
            }
        }
        for(auto v: mani.vertices())
            mani.pos(v) = xform.inverse(mani.pos(v));
    }
    
}
