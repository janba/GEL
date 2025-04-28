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
#include <GEL/HMesh/dual.h>

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

namespace  HMesh {

    HMesh::Manifold volume_polygonize(const XForm& xform, const Geometry::RGrid<float>& grid,
                                      float tau, bool make_triangles, bool high_is_inside,
                                      bool dual_connectivity)
    {
        Manifold mani;
        
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
        
        vector<Vec3d> edge_intersections;
        for(Vec3i pi: Range3D(grid.get_dims()))
            if(is_inside(pi)) {
                Vec3d p(pi);
                for (int nbr_idx = 0; nbr_idx < 6 ; ++ nbr_idx) {
                    Vec3i pni = pi + N6i[nbr_idx];
                    if(is_outside(pni)) {
                        float t = 0;
                        if (grid.in_domain(pni)) {
                            float va = grid[pi];
                            float vb = grid[pni];
                            t = (tau-va)/(vb - va);
                        }
                        edge_intersections.push_back(p * (1-t) + Vec3d(pni) * t);
                        vector<Vec3d> quad_vertices(4);
                        for(int n=0;n<4;++n)
                            quad_vertices[n] = p + hex_faces[nbr_idx][3-n];
                        mani.add_face(quad_vertices);
                        
                    }
                }
            }
                
        stitch_mesh(mani, 1e-10);
        
        if (dual_connectivity) {
            for(auto v: mani.vertices()) {
                Vec3d p(0);
                int cnt = 0;
                for (auto f: mani.incident_faces(v)) 
                    if (mani.in_use(f)) {
                        auto id = f.get_index();
                        p += edge_intersections[id];
                        ++cnt;
                }
                mani.pos(v) = p/cnt;
            }
        }
        else {
            Manifold mc_mesh;
            VertexAttributeVector<int> cluster_id(-1);
            for (auto v: mani.vertices()) {
                vector<Vec3d> points;
                vector<int> ids;
                ids.reserve(20);
                for (auto f: mani.incident_faces(v))
                    if (mani.in_use(f)) {
                        auto id = f.get_index();
                        points.push_back(edge_intersections[id]);
                        ids.push_back(id);
                    }
                if(points.size()>2) {
                    FaceID f_mc = mc_mesh.add_face(points);
                    int i = 0;
                    for (auto vv: mc_mesh.incident_vertices(f_mc))
                        cluster_id[vv] = ids[i++];
                }
            }
            stitch_mesh(mc_mesh, cluster_id);
            mani = mc_mesh;
        }
        
        for(auto v: mani.vertices())
            mani.pos(v) = xform.inverse(mani.pos(v));
        
        if(make_triangles) {
            triangulate(mani, SHORTEST_EDGE);
            remove_needles(mani, 0.1, true);
        }
        
        mani.cleanup();
        
        return mani;
    }

}
