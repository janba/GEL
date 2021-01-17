//
//  polygonize.cpp
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

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

    
    void polygonize(const XForm& xform, const RGrid<float>& grid, std::vector<CGLA::Vec3d>& quad_vertices,
                    float tau, bool high_is_inside)
    {
        auto is_inside = [&](const Vec3i& pi) {
            return high_is_inside == (grid[pi] > tau);
        };
        auto is_outside = [&](const Vec3i& pi) {
            return !grid.in_domain(pi) || (high_is_inside == (grid[pi] <= tau));
        };
        
        quad_vertices.clear();
        for(int i=0;i<xform.get_dims()[0];++i)
            for(int j=0;j<xform.get_dims()[1];++j)
                for(int k=0;k<xform.get_dims()[2];++k)
                {
                    Vec3i pi(i,j,k);
                    Vec3d p(pi);
                    if(is_inside(pi))
                    {
                        for (int hf_idx = 0; hf_idx < 6 ; ++ hf_idx) {
                            Vec3i pni = pi + N6i[hf_idx];
                            if(is_outside(pni)) {
                                Vec3d pn = p + 0.5 * N6d[hf_idx];
                                if(high_is_inside)
                                    for(int n=0;n<4;++n)
                                        quad_vertices.push_back(xform.inverse(pn + hex_faces[hf_idx][3-n]));
                                else
                                    for(int n=0;n<4;++n)
                                        quad_vertices.push_back(xform.inverse(pn + hex_faces[hf_idx][n]));
                            }
                        }
                    }
                }
        
    }


    bool towards_intersection(const XForm& xform, const Geometry::RGrid<float>& grid, Vec3d& p, const Vec3d& p1, double tau) {
        auto v0 = interpolate(grid,xform.apply(p));
        auto v1 = interpolate(grid,xform.apply(p1));
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
        const double delta = (sqrt(3.0)/2.0) * xform.inv_scale();

        mani.clear();
        vector<Vec3d> quad_vertices;
        polygonize(xform, grid, quad_vertices, tau, high_is_inside);
        vector<int> indices;
        vector<int> faces(quad_vertices.size()/4,4);
        for(int i=0;i<quad_vertices.size();++i)
            indices.push_back(i);
        build(mani, quad_vertices.size(),
                   quad_vertices[0].get(),
                   faces.size(),
                   &faces[0],
                   &indices[0]);
        
        stitch_more(mani, 0.001 * delta);
        
        if(make_triangles)
            triangulate(mani);
        
        mani.cleanup();
        
        laplacian_smooth(mani, .25, 3);
        for(auto v: mani.vertices()) {
            Vec3d& p = mani.pos(v);
            const Vec3d n = normal(mani,v);
            const Vec3d p1 = p + delta * n;
            if(!towards_intersection(xform, grid, p, p1, tau)) {
                const Vec3d p1 = p - delta * n;
                towards_intersection(xform, grid, p, p1, tau);
            }
        }
    }
    
}
