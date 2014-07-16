//
//  polygonize.cpp
//  PointReconstruction
//
//  Created by J. Andreas Bærentzen on 16/03/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#include "../Geometry/Implicit.h"
#include "polygonize.h"
#include "smooth.h"
#include "cleanup.h"
#include "triangulate.h"

using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace HMesh;

namespace {
    Vec3d xpf[] = {Vec3d(0,-0.5,-0.5),Vec3d(0,0.5,-0.5),Vec3d(0,0.5,0.5),Vec3d(0,-0.5,0.5)};
    Vec3d xmf[] = {Vec3d(0, 0.5,-0.5),Vec3d(0,-0.5,-0.5),Vec3d(0,-0.5,0.5),Vec3d(0,0.5,0.5)};
    
    Vec3d ypf[] = {Vec3d( 0.5,0, -0.5),Vec3d(-0.5,0, -0.5),Vec3d(-0.5,0, 0.5),Vec3d(0.5,0, 0.5)};
    Vec3d ymf[] = {Vec3d(-0.5,0, -0.5),Vec3d(0.5,0, -0.5),Vec3d(0.5,0, 0.5),Vec3d(-0.5,0, 0.5)};
    
    Vec3d zpf[] = {Vec3d(-0.5,-0.5,0),Vec3d(0.5,-0.5,0),Vec3d(0.5,0.5,0),Vec3d(-0.5,0.5,0)};
    Vec3d zmf[] = {Vec3d( 0.5,-0.5,0),Vec3d(-0.5,-0.5,0),Vec3d(-0.5,0.5,0),Vec3d(0.5,0.5,0)};
}

namespace HMesh
{
    
    void polygonize(const XForm& xform, const RGrid<float>& grid, std::vector<CGLA::Vec3d>& quad_vertices, float tau)
    {
        quad_vertices.clear();
        for(int i=0;i<xform.get_dims()[0];++i)
            for(int j=0;j<xform.get_dims()[1];++j)
                for(int k=0;k<xform.get_dims()[2];++k)
                {
                    Vec3i vox(i,j,k);
                    if(grid[vox] <= tau)
                    {
                        if(grid.in_domain(Vec3i(i+1,j,k)) && grid[Vec3i(i+1,j,k)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i+0.5, j,k) + xpf[n]));
                        if(grid.in_domain(Vec3i(i-1,j,k)) && grid[Vec3i(i-1,j,k)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i-0.5, j,k) + xmf[n]));
                        if(grid.in_domain(Vec3i(i,j+1,k)) && grid[Vec3i(i,j+1,k)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i, j+0.5,k) + ypf[n]));
                        if(grid.in_domain(Vec3i(i,j-1,k)) && grid[Vec3i(i,j-1,k)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i, j-0.5,k) + ymf[n]));
                        if(grid.in_domain(Vec3i(i,j,k+1)) && grid[Vec3i(i,j,k+1)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i, j,k+0.5) + zpf[n]));
                        if(grid.in_domain(Vec3i(i,j,k-1)) && grid[Vec3i(i,j,k-1)] > tau)
                            for(int n=0;n<4;++n)
                                quad_vertices.push_back(xform.inverse(Vec3d(i, j,k-0.5) + zmf[n]));
                        
                    }
                }
        
    }
    
    
    void volume_polygonize(const XForm& xform, const Geometry::RGrid<float>& grid,
                           HMesh::Manifold& mani, float tau)
    {
        mani.clear();
        vector<Vec3d> quad_vertices;
        polygonize(xform, grid, quad_vertices, tau);
        vector<int> indices;
        vector<int> faces(quad_vertices.size()/4,4);
        for(int i=0;i<quad_vertices.size();++i)
            indices.push_back(i);
        mani.build(quad_vertices.size(),
                   quad_vertices[0].get(),
                   faces.size(),
                   &faces[0],
                   &indices[0]);
        
        stitch_more(mani, 1e-5);
        shortest_edge_triangulate(mani);
        
        float avg_edge_len=0;
        for(HalfEdgeIDIterator h = mani.halfedges_begin(); h != mani.halfedges_end();++h)
            avg_edge_len += length(mani, *h);
        avg_edge_len /= mani.no_halfedges();
        VolumetricImplicit imp(xform, grid);
        for(int iter=0;iter<4;++iter)
        {
            TAL_smoothing(mani, .25, 1);
            for(VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
                imp.push_to_surface(mani.pos(*vid),0,avg_edge_len*0.5);
        }
        mani.cleanup();
        cout << "Produced" << mani.no_faces() << " faces " << endl;
    }
    
}
