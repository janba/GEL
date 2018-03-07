//
//  hmesh_functions.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef hmesh_functions_hpp
#define hmesh_functions_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <GEL/HMesh/HMesh.h>

extern "C" {
    DLLEXPORT void stitch_mesh(HMesh::Manifold* m_ptr, double rad);
    
    DLLEXPORT bool valid(const HMesh::Manifold* m_ptr);
    DLLEXPORT bool closed(const HMesh::Manifold* m_ptr);
    
    DLLEXPORT void bbox(const HMesh::Manifold* m_ptr, CGLA::Vec3d* pmin, CGLA::Vec3d* pmax);
    DLLEXPORT void bsphere(const HMesh::Manifold* m_ptr, CGLA::Vec3d* c, double* r);

    DLLEXPORT bool obj_load(char*, HMesh::Manifold*);
    DLLEXPORT bool off_load(char*, HMesh::Manifold* m_ptr);
    DLLEXPORT bool ply_load(char*, HMesh::Manifold* m_ptr);
    DLLEXPORT bool x3d_load(char*, HMesh::Manifold* m_ptr);
    
    DLLEXPORT bool obj_save(char*, HMesh::Manifold* m_ptr);
    DLLEXPORT bool off_save(char*, HMesh::Manifold* m_ptr);
    DLLEXPORT bool x3d_save(char*, HMesh::Manifold* m_ptr);

        
    DLLEXPORT void remove_caps(HMesh::Manifold* m_ptr, float thresh);
    
    DLLEXPORT void remove_needles(HMesh::Manifold* m_ptr, float thresh, bool averagePositions = false);

    DLLEXPORT void close_holes(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void flip_orientation(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void minimize_curvature(HMesh::Manifold* m_ptr, bool anneal=false);
    
    DLLEXPORT void minimize_dihedral_angle(HMesh::Manifold* m_ptr, int max_iter=10000, bool anneal=false, bool alpha=false, double gamma=4.0);
    
    DLLEXPORT void maximize_min_angle(HMesh::Manifold* m_ptr, float thresh, bool anneal=false);
    
    DLLEXPORT void optimize_valency(HMesh::Manifold* m_ptr, bool anneal=false);
    
    DLLEXPORT void randomize_mesh(HMesh::Manifold* m_ptr, int max_iter);
    
    DLLEXPORT void quadric_simplify(HMesh::Manifold* m_ptr, double keep_fraction, double singular_thresh = 0.0001, bool choose_optimal_positions = true);

    DLLEXPORT float average_edge_length(const HMesh::Manifold* m_ptr);
    
    DLLEXPORT float median_edge_length(const HMesh::Manifold* m_ptr);
    
    DLLEXPORT int refine_edges(HMesh::Manifold* m_ptr, float t);

    DLLEXPORT void cc_split(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void loop_split(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void root3_subdivide(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void rootCC_subdivide(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void butterfly_subdivide(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void cc_smooth(HMesh::Manifold* m_ptr);
    
    DLLEXPORT void loop_smooth(HMesh::Manifold* m_ptr);

    DLLEXPORT void shortest_edge_triangulate(HMesh::Manifold* m_ptr);

}

#endif /* hmesh_functions_hpp */
