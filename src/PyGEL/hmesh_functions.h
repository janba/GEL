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

#include <stdbool.h>

typedef  char* Manifold_ptr;

#ifdef __cplusplus
extern "C" {
#endif
    DLLEXPORT int stitch_mesh(Manifold_ptr m_ptr, double rad);
    
    DLLEXPORT bool valid(const Manifold_ptr m_ptr);
    DLLEXPORT bool closed(const Manifold_ptr m_ptr);
    
    DLLEXPORT void bbox(const Manifold_ptr m_ptr, double* pmin, double* pmax);
    DLLEXPORT void bsphere(const Manifold_ptr m_ptr, double* c, double* r);

    DLLEXPORT bool obj_load(const char*, Manifold_ptr);
    DLLEXPORT bool off_load(const char*, Manifold_ptr m_ptr);
    DLLEXPORT bool ply_load(const char*, Manifold_ptr m_ptr);
    DLLEXPORT bool x3d_load(const char*, Manifold_ptr m_ptr);
    
    DLLEXPORT bool obj_save(const char*, Manifold_ptr m_ptr);
    DLLEXPORT bool off_save(const char*, Manifold_ptr m_ptr);
    DLLEXPORT bool x3d_save(const char*, Manifold_ptr m_ptr);

        
    DLLEXPORT void remove_caps(Manifold_ptr m_ptr, float thresh);
    
    DLLEXPORT void remove_needles(Manifold_ptr m_ptr, float thresh, bool averagePositions);

    DLLEXPORT void close_holes(Manifold_ptr m_ptr, int max_size);
    
    DLLEXPORT void flip_orientation(Manifold_ptr m_ptr);
    
    DLLEXPORT void merge_coincident_boundary_vertices(Manifold_ptr m_ptr, double rad);
    
    DLLEXPORT void minimize_curvature(Manifold_ptr m_ptr, bool anneal);
    
    DLLEXPORT void minimize_dihedral_angle(Manifold_ptr m_ptr, int max_iter, bool anneal, bool alpha, double gamma);
    
    DLLEXPORT void maximize_min_angle(Manifold_ptr m_ptr, float thresh, bool anneal);
    
    DLLEXPORT void optimize_valency(Manifold_ptr m_ptr, bool anneal);
    
    DLLEXPORT void randomize_mesh(Manifold_ptr m_ptr, int max_iter);
    
    DLLEXPORT void quadric_simplify(Manifold_ptr m_ptr, double keep_fraction, double singular_thresh, bool choose_optimal_positions);

    DLLEXPORT float average_edge_length(const Manifold_ptr m_ptr);
    
    DLLEXPORT float median_edge_length(const Manifold_ptr m_ptr);
    
    DLLEXPORT int refine_edges(Manifold_ptr m_ptr, float t);

    DLLEXPORT void cc_split(Manifold_ptr m_ptr);
    
    DLLEXPORT void loop_split(Manifold_ptr m_ptr);
    
    DLLEXPORT void root3_subdivide(Manifold_ptr m_ptr);
    
    DLLEXPORT void rootCC_subdivide(Manifold_ptr m_ptr);
    
    DLLEXPORT void butterfly_subdivide(Manifold_ptr m_ptr);
    
    DLLEXPORT void cc_smooth(Manifold_ptr m_ptr);
    
    DLLEXPORT void loop_smooth(Manifold_ptr m_ptr);

    DLLEXPORT void shortest_edge_triangulate(Manifold_ptr m_ptr);

    DLLEXPORT void ear_clip_triangulate(Manifold_ptr m_ptr);

#ifdef __cplusplus
}
#endif

#endif /* hmesh_functions_hpp */
