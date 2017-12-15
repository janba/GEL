//
//  hmesh_functions.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef hmesh_functions_hpp
#define hmesh_functions_hpp

#include <GEL/HMesh/HMesh.h>

extern "C" {
    void stitch_mesh(HMesh::Manifold* m_ptr, double rad);
    
    bool valid(const HMesh::Manifold* m_ptr);
    bool closed(const HMesh::Manifold* m_ptr);
    
    void bbox(const HMesh::Manifold* m_ptr, CGLA::Vec3d* pmin, CGLA::Vec3d* pmax);
    void bsphere(const HMesh::Manifold* m_ptr, CGLA::Vec3d* c, double* r);

    bool obj_load(char*, HMesh::Manifold*);
    bool off_load(char*, HMesh::Manifold* m_ptr);
    bool ply_load(char*, HMesh::Manifold* m_ptr);
    bool x3d_load(char*, HMesh::Manifold* m_ptr);
        
    void remove_caps(HMesh::Manifold* m_ptr, float thresh);
    
    void remove_needles(HMesh::Manifold* m_ptr, float thresh, bool averagePositions = false);

    void close_holes(HMesh::Manifold* m_ptr);
    
    void flip_orientation(HMesh::Manifold* m_ptr);
    
    void minimize_curvature(HMesh::Manifold* m_ptr, bool anneal=false);
    
    void maximize_min_angle(HMesh::Manifold* m_ptr, float thresh, bool anneal=false);
    
    void optimize_valency(HMesh::Manifold* m_ptr, bool anneal=false);
    
    void randomize_mesh(HMesh::Manifold* m_ptr, int max_iter);
    
    void quadric_simplify(HMesh::Manifold* m_ptr, double keep_fraction, double singular_thresh = 0.0001, bool choose_optimal_positions = true);

    float average_edge_length(const HMesh::Manifold* m_ptr);
    
    float median_edge_length(const HMesh::Manifold* m_ptr);
    
    int refine_edges(HMesh::Manifold* m_ptr, float t);

    void cc_split(HMesh::Manifold* m_ptr);
    
    void loop_split(HMesh::Manifold* m_ptr);
    
    void root3_subdivide(HMesh::Manifold* m_ptr);
    
    void rootCC_subdivide(HMesh::Manifold* m_ptr);
    
    void butterfly_subdivide(HMesh::Manifold* m_ptr);
    
    void cc_smooth(HMesh::Manifold* m_ptr);
    
    void loop_smooth(HMesh::Manifold* m_ptr);

    void shortest_edge_triangulate(HMesh::Manifold* m_ptr);

}

#endif /* hmesh_functions_hpp */
