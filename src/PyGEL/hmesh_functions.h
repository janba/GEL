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

    
    int stitch_mesh(HMesh::Manifold* m_ptr, double rad);
    bool valid(const HMesh::Manifold* m_ptr);
    bool closed(const HMesh::Manifold* m_ptr);

    std::pair<std::vector<double>, std::vector<double>> bbox(const HMesh::Manifold* m_ptr);
    std::pair<std::vector<double>, double> bsphere(const HMesh::Manifold* m_ptr);

    bool load(const std::string& filename, HMesh::Manifold* m_ptr);
    bool obj_load(const std::string& filename, HMesh::Manifold* m_ptr);
    bool off_load(const std::string& filename, HMesh::Manifold* m_ptr);
    bool ply_load(const std::string& filename, HMesh::Manifold* m_ptr);
    bool x3d_load(const std::string& filename, HMesh::Manifold* m_ptr);

    bool obj_save(const std::string& filename, HMesh::Manifold* m_ptr);
    bool off_save(const std::string& filename, HMesh::Manifold* m_ptr);
    bool x3d_save(const std::string& filename, HMesh::Manifold* m_ptr);

    void remove_caps(HMesh::Manifold* m_ptr, float thresh);
    void remove_needles(HMesh::Manifold* m_ptr, float thresh, bool averagePositions);
    void close_holes(HMesh::Manifold* m_ptr, int max_size);
    void flip_orientation(HMesh::Manifold* m_ptr);
    void merge_coincident_boundary_vertices(HMesh::Manifold* m_ptr, double rad);

    void minimize_curvature(HMesh::Manifold* m_ptr, bool anneal);
    void minimize_dihedral_angle(HMesh::Manifold* m_ptr, int max_iter, bool anneal, bool alpha, double gamma);
    void maximize_min_angle(HMesh::Manifold* m_ptr, float thresh, bool anneal);
    void optimize_valency(HMesh::Manifold* m_ptr, bool anneal);
    void randomize_mesh(HMesh::Manifold* m_ptr, int max_iter);

    void quadric_simplify(HMesh::Manifold* m_ptr, double keep_fraction, double singular_thresh, double error_thresh);
    float average_edge_length(const HMesh::Manifold* m_ptr);
    float median_edge_length(const HMesh::Manifold* m_ptr);
    int refine_edges(HMesh::Manifold* m_ptr, float t);

    void cc_split(HMesh::Manifold* m_ptr);
    void loop_split(HMesh::Manifold* m_ptr);
    void root3_subdivide(HMesh::Manifold* m_ptr);
    void rootCC_subdivide(HMesh::Manifold* m_ptr);
    void butterfly_subdivide(HMesh::Manifold* m_ptr);

    void cc_smooth(HMesh::Manifold* m_ptr);
    void volume_preserving_cc_smooth(HMesh::Manifold* m_ptr, int iter);
    void regularize_quads(HMesh::Manifold* m_ptr, float weight, float shrink, int iter);
    void loop_smooth(HMesh::Manifold* m_ptr);
    void taubin_smooth(HMesh::Manifold* m_ptr, int iter);
    void laplacian_smooth(HMesh::Manifold* m_ptr, float weight, int iter);
    void anisotropic_smooth(HMesh::Manifold* m_ptr, float sharpness, int iter);

    void volumetric_isocontour(HMesh::Manifold* m_ptr, int x_dim, int y_dim, int z_dim, const std::vector<float>& data,
                               const std::vector<double>& pmin, const std::vector<double>& pmax,
                               float tau, bool make_triangles, bool high_is_inside, bool dual_connectivity);

    void shortest_edge_triangulate(HMesh::Manifold* m_ptr);
    void ear_clip_triangulate(HMesh::Manifold* m_ptr);

    void graph_to_feq(Geometry::AMGraph3D* g_ptr, HMesh::Manifold* m_ptr, const std::vector<double>& node_radii, bool symmetrize, bool use_graph_radii);
    void non_rigid_registration(HMesh::Manifold* m_ptr, HMesh::Manifold* m_ref_ptr);

    void rsr_recon(HMesh::Manifold* m_ptr, const std::vector<double>& verts, const std::vector<double>& normals, 
                   int v_num, int n_num, bool isEuclidean = false, int genus = 0, int k = 70, int r = 20, int theta = 60, int n = 50);

    void extrude_faces(HMesh::Manifold* m_ptr, const std::vector<int>& faces, std::vector<size_t>& fidx_ptr);
    void kill_face_loop(HMesh::Manifold* m_ptr);
    void kill_degenerate_face_loops(HMesh::Manifold* m_ptr, double thresh);
    void stable_marriage_registration(HMesh::Manifold* m_ptr, HMesh::Manifold* m_ref_ptr);

    std::vector<HMesh::Manifold*>* connected_components(HMesh::Manifold* m_ptr);
    size_t mesh_vec_size(const std::vector<HMesh::Manifold*>* mv_ptr);
    HMesh::Manifold* mesh_vec_get(const std::vector<HMesh::Manifold*>* mv_ptr, size_t i);
    void mesh_vec_del(std::vector<HMesh::Manifold*>* mv_ptr);
    
    int count_boundary_curves(HMesh::Manifold* m_ptr);
    // void create_LBO(HMesh::Manifold* m_ptr, std::vector<size_t>& L_i, std::vector<size_t>& L_j, std::vector<double>& L,
    //                 std::vector<size_t>& M_i, std::vector<size_t>& M_j, std::vector<double>& M);
}

#endif /* hmesh_functions_hpp */
