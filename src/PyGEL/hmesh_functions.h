//
//  hmesh_functions.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef hmesh_functions_hpp
#define hmesh_functions_hpp

#include <vector>
#include <string>
#include "Vec3dVector.h"

#include "Manifold.h"
#include "Graph.h"

namespace PyGEL {

    
    int stitch_mesh(Manifold* m_ptr, double rad);
    bool valid(const Manifold* m_ptr);
    bool closed(const Manifold* m_ptr);

    std::pair<std::vector<double>, std::vector<double>> bbox(const Manifold* m_ptr);
    std::pair<std::vector<double>, double> bsphere(const Manifold* m_ptr);

    bool load(const std::string& filename, Manifold* m_ptr);
    bool obj_load(const std::string& filename, Manifold* m_ptr);
    bool off_load(const std::string& filename, Manifold* m_ptr);
    bool ply_load(const std::string& filename, Manifold* m_ptr);
    bool x3d_load(const std::string& filename, Manifold* m_ptr);

    bool obj_save(const std::string& filename, Manifold* m_ptr);
    bool off_save(const std::string& filename, Manifold* m_ptr);
    bool x3d_save(const std::string& filename, Manifold* m_ptr);

    void remove_caps(Manifold* m_ptr, float thresh);
    void remove_needles(Manifold* m_ptr, float thresh, bool averagePositions);
    void close_holes(Manifold* m_ptr, int max_size);
    void flip_orientation(Manifold* m_ptr);
    void merge_coincident_boundary_vertices(Manifold* m_ptr, double rad);

    void minimize_curvature(Manifold* m_ptr, bool anneal);
    void minimize_dihedral_angle(Manifold* m_ptr, int max_iter, bool anneal, bool alpha, double gamma);
    void maximize_min_angle(Manifold* m_ptr, float thresh, bool anneal);
    void optimize_valency(Manifold* m_ptr, bool anneal);
    void randomize_mesh(Manifold* m_ptr, int max_iter);

    void quadric_simplify(Manifold* m_ptr, double keep_fraction, double singular_thresh, double error_thresh);
    float average_edge_length(const Manifold* m_ptr);
    float median_edge_length(const Manifold* m_ptr);
    int refine_edges(Manifold* m_ptr, float t);

    void cc_split(Manifold* m_ptr);
    void loop_split(Manifold* m_ptr);
    void root3_subdivide(Manifold* m_ptr);
    void rootCC_subdivide(Manifold* m_ptr);
    void butterfly_subdivide(Manifold* m_ptr);

    void cc_smooth(Manifold* m_ptr);
    void volume_preserving_cc_smooth(Manifold* m_ptr, int iter);
    void regularize_quads(Manifold* m_ptr, float weight, float shrink, int iter);
    void loop_smooth(Manifold* m_ptr);
    void taubin_smooth(Manifold* m_ptr, int iter);
    void laplacian_smooth(Manifold* m_ptr, float weight, int iter);
    void anisotropic_smooth(Manifold* m_ptr, float sharpness, int iter);

    void volumetric_isocontour(Manifold* m_ptr, int x_dim, int y_dim, int z_dim, const std::vector<float>& data,
                               const std::vector<double>& pmin, const std::vector<double>& pmax,
                               float tau, bool make_triangles, bool high_is_inside, bool dual_connectivity);

    void shortest_edge_triangulate(Manifold* m_ptr);
    void ear_clip_triangulate(Manifold* m_ptr);

    void graph_to_feq(AMGraph3D* g_ptr, Manifold* m_ptr, const std::vector<double>& node_radii, bool symmetrize, bool use_graph_radii);
    void non_rigid_registration(Manifold* m_ptr, Manifold* m_ref_ptr);

    void rsr_recon(Manifold* m_ptr, const std::vector<double>& verts, const std::vector<double>& normals, 
                   int v_num, int n_num, bool isEuclidean = false, int genus = 0, int k = 70, int r = 20, int theta = 60, int n = 50);

    void extrude_faces(Manifold* m_ptr, const std::vector<int>& faces, std::vector<size_t>& fidx_ptr);
    void kill_face_loop(Manifold* m_ptr);
    void kill_degenerate_face_loops(Manifold* m_ptr, double thresh);
    void stable_marriage_registration(Manifold* m_ptr, Manifold* m_ref_ptr);

    std::vector<Manifold*>* connected_components(Manifold* m_ptr);
    size_t mesh_vec_size(const std::vector<Manifold*>* mv_ptr);
    Manifold* mesh_vec_get(const std::vector<Manifold*>* mv_ptr, size_t i);
    void mesh_vec_del(std::vector<Manifold*>* mv_ptr);
    
    int count_boundary_curves(Manifold* m_ptr);
    // void create_LBO(Manifold* m_ptr, std::vector<size_t>& L_i, std::vector<size_t>& L_j, std::vector<double>& L,
    //                 std::vector<size_t>& M_i, std::vector<size_t>& M_j, std::vector<double>& M);
}

#endif /* hmesh_functions_hpp */
