//
//  Delaunay_triangulate.cpp
//  GEL
//
//  Created by Jakob Andreas Bærentzen on 13/09/2018.
//  Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
//

#include "../Geometry/jrs_triangle.h"
#include "Delaunay_triangulate.h"
#include "HMesh.h"
#include "../CGLA/Vec2d.h"

using namespace CGLA;
using namespace std;

namespace HMesh {
    
    HMesh::Manifold Delaunay_triangulate(const std::vector<CGLA::Vec3d>& pts3d, const CGLA::Vec3d& X_axis, const CGLA::Vec3d& Y_axis) {
        vector<Vec2d> pts2d(pts3d.size());
        for (int i=0;i<pts3d.size();++i)
            pts2d[i] = Vec2d(dot(pts3d[i], X_axis), dot(pts3d[i], Y_axis));
        
        // The code below builds the triangulation
        triangulateio pts_in, tri_out;
        pts_in.numberofpoints = pts2d.size();
        pts_in.pointlist = static_cast<double*>(&(pts2d[0][0]));
        pts_in.numberofpointattributes = 0;
        pts_in.pointmarkerlist = 0;
        pts_in.segmentmarkerlist = 0;
        pts_in.numberofholes = 0;
        pts_in.numberofregions = 0;
        
        tri_out.pointlist = 0;
        tri_out.pointmarkerlist = 0;
        tri_out.pointattributelist = 0;
        tri_out.numberofpointattributes = 0;
        tri_out.segmentlist = 0;
        tri_out.trianglelist = 0;
        
        // Call Triangle with arguments that specify: (z)ero is firt index, (p)slg triangulation.
        // no (B)oundary markers, no (S)teiner points not absolutely needed. Operate (Q)ietly.
        string triangulate_cmd_str("zBSQ");
        triangulate("zBSQ", &pts_in, &tri_out, 0);
        
        if(tri_out.numberofpoints > pts_in.numberofpoints){
            cout << "Steiner points were created, any incident triangles will be removed..." << endl;
        }
        
        // Now, go through all triangles created and add them to the mesh. If a vertex of a triangle
        // is detected as belonging to the ROI boundary or a forced Steiner point, it is removed.
        // If a triangle belongs to the exterior of the boundary, we also remove it. This test is tricky,
        // but we can use the fact that edges are oriented. If all three vertices belong to the boundary,
        // the outside triangles have two edges with same orientation as the boundary chain orientation
        // whereas the interior triangles have only one.
        
        Manifold m_new;
        VertexAttributeVector<int> vid_to_index;
        
        for(int l=0;l<tri_out.numberoftriangles;++l) {
            int i = tri_out.trianglelist[3*l+0];
            int j = tri_out.trianglelist[3*l+1];
            int k = tri_out.trianglelist[3*l+2];
            
            if(!(i<pts2d.size() && j<pts2d.size() && k<pts2d.size())) {
//                cout << "Removing face incident on Steiner " << endl;
                continue;
            }
            
            FaceID f = m_new.add_face({pts3d[i], pts3d[j], pts3d[k]});
            Walker w = m_new.walker(f);
            vid_to_index[w.vertex()] = i; w = w.next();
            vid_to_index[w.vertex()] = j; w = w.next();
            vid_to_index[w.vertex()] = k;
        }
        
        // Release memory used by triangle.
        trifree(tri_out.trianglelist);
        trifree(tri_out.pointlist);
        trifree(tri_out.pointattributelist);
        trifree(tri_out.segmentlist);
        trifree(tri_out.pointmarkerlist);
        
        // Stitch the added faces. the mapping to global index is used to sort out connectivity.
        stitch_mesh(m_new,vid_to_index);
        
        return m_new;
    }
}
