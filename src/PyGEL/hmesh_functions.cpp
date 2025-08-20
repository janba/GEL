//
//  hmesh_functions.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 04/11/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "hmesh_functions.h"
#include <string>
#include <GEL/HMesh/HMesh.h>
#include <GEL/HMesh/face_loop.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_skeletonize.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/HMesh/RsR.h>
#include "Graph.h"
#include "Manifold.h"



using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Geometry;

namespace PyGEL {

    int stitch_mesh(Manifold* m_ptr, double rad) {
        return stitch_mesh(*m_ptr, rad);
    }

    bool valid(const Manifold* m_ptr) {
        return valid(*m_ptr);
    }

    bool closed(const Manifold* m_ptr) {
        return closed(*m_ptr);
    }

    bool load(const std::string& filename, Manifold* m_ptr) {
        return load(filename, *m_ptr);
    }
    bool obj_load(const std::string& filename, Manifold* m_ptr) {
        return obj_load(filename, *m_ptr);
    }
    bool off_load(const std::string& filename, Manifold* m_ptr) {
        return off_load(filename, *m_ptr);
    }
    bool ply_load(const std::string& filename, Manifold* m_ptr) {
        return ply_load(filename, *m_ptr);
    }
    bool x3d_load(const std::string& filename, Manifold* m_ptr) {
        return x3d_load(filename, *m_ptr);
    }

    bool obj_save(const std::string& filename, Manifold* m_ptr) {
        return obj_save(filename, *m_ptr);
    }
    bool off_save(const std::string& filename, Manifold* m_ptr) {
        return off_save(filename, *m_ptr);
    }
    bool x3d_save(const std::string& filename, Manifold* m_ptr) {
        return x3d_save(filename, *m_ptr);
    }

    std::pair<std::vector<double>, std::vector<double>> bbox(const Manifold* m_ptr) {
        Vec3d pmin, pmax;
        bbox(*m_ptr, pmin, pmax);
        std::vector<double> vmin{pmin[0], pmin[1], pmin[2]};
        std::vector<double> vmax{pmax[0], pmax[1], pmax[2]};
        return {vmin, vmax};
    }

    std::pair<std::vector<double>, double> bsphere(const Manifold* m_ptr) {
        Vec3d c; float r = 0.0f;
        bsphere(*m_ptr, c, r);
        std::vector<double> vc{c[0], c[1], c[2]};
        return {vc, static_cast<double>(r)};
    }


    void remove_caps(Manifold* m_ptr, float thresh) {
        remove_caps(*m_ptr, thresh);
    }

    void remove_needles(Manifold* m_ptr, float thresh, bool averagePositions) {
        remove_needles(*m_ptr, thresh, averagePositions);
    }

    void close_holes(Manifold* m_ptr, int max_size) {
        close_holes(*m_ptr, max_size);
    }

    void flip_orientation(Manifold* m_ptr) {
        flip_orientation(*m_ptr);
    }

    void merge_coincident_boundary_vertices(Manifold* m_ptr, double rad) {
        merge_coincident_boundary_vertices(*m_ptr, rad);
    }

    void minimize_curvature(Manifold* m_ptr, bool anneal) {
        minimize_curvature(*m_ptr, anneal);
    }

    void minimize_dihedral_angle(Manifold* m_ptr, int max_iter, bool anneal, bool alpha, double gamma) {
        minimize_dihedral_angle(*m_ptr, max_iter, anneal, alpha, gamma);
    }

    void maximize_min_angle(Manifold* m_ptr, float thresh, bool anneal) {
        maximize_min_angle(*m_ptr, thresh, anneal);
    }

    void optimize_valency(Manifold* m_ptr, bool anneal) {
        optimize_valency(*m_ptr, anneal);
    }

    void randomize_mesh(Manifold* m_ptr, int max_iter) {
        randomize_mesh(*m_ptr, max_iter);
    }

    void quadric_simplify(Manifold* m_ptr, double keep_fraction,
                        double singular_thresh,
                        double error_thresh) {
        quadric_simplify(*m_ptr, keep_fraction, singular_thresh, error_thresh);
    }

    float average_edge_length(const Manifold* m_ptr) {
        return average_edge_length(*m_ptr);
    }

    float median_edge_length(const Manifold* m_ptr) {
        return median_edge_length(*m_ptr);
    }

    int refine_edges(Manifold* m_ptr, float t) {
        return refine_edges(*m_ptr, t);
    }

    void cc_split(Manifold* m_ptr) {
        cc_split(*m_ptr);
    }

    void loop_split(Manifold* m_ptr) {
        loop_split(*m_ptr);
    }

    void root3_subdivide(Manifold* m_ptr) {
        root3_subdivide(*m_ptr);
    }

    void rootCC_subdivide(Manifold* m_ptr) {
        rootCC_subdivide(*m_ptr);
    }

    void butterfly_subdivide(Manifold* m_ptr) {
        butterfly_subdivide(*m_ptr);
    }

    void cc_smooth(Manifold* m_ptr) {
        cc_smooth(*m_ptr);
    }

    void volume_preserving_cc_smooth(Manifold* m_ptr, int iter) {
        volume_preserving_cc_smooth(*m_ptr, iter);
    }

    void regularize_quads(Manifold* m_ptr, float weight, float shrink, int iter) {
        regularize_quads(*m_ptr, weight, shrink, iter);
    }

    void loop_smooth(Manifold* m_ptr) {
        loop_smooth(*m_ptr);
    }

    void shortest_edge_triangulate(Manifold* m_ptr) {
        triangulate(*m_ptr, SHORTEST_EDGE);
    }

    void ear_clip_triangulate(Manifold* m_ptr) {
        triangulate(*m_ptr, CLIP_EAR);
    }

    void taubin_smooth(Manifold* m_ptr, int iter) {
        taubin_smooth(*m_ptr, iter);
    }

    void laplacian_smooth(Manifold* m_ptr, float weight, int iter) {
        laplacian_smooth(*m_ptr, weight, iter);
    }

    void anisotropic_smooth(Manifold* m_ptr, float sharpness, int iter) {
        anisotropic_smooth(*m_ptr, iter, sharpness);
    }

    void volumetric_isocontour(Manifold* m_ptr, int x_dim, int y_dim, int z_dim,
                            const std::vector<float>& data,
                            const std::vector<double>& pmin_v,
                            const std::vector<double>& pmax_v,
                            float tau, bool make_triangles,
                            bool high_is_inside,
                            bool dual_connectivity) {
        Vec3i dims(x_dim, y_dim, z_dim);
        Vec3d pmin(pmin_v.size() > 0 ? pmin_v[0] : 0.0,
                pmin_v.size() > 1 ? pmin_v[1] : 0.0,
                pmin_v.size() > 2 ? pmin_v[2] : 0.0);
        Vec3d pmax(pmax_v.size() > 0 ? pmax_v[0] : 0.0,
                pmax_v.size() > 1 ? pmax_v[1] : 0.0,
                pmax_v.size() > 2 ? pmax_v[2] : 0.0);
        XForm xform(pmin, pmax, dims, 0.0);
        RGrid<float> grid(dims);
        // Expect data.size() == grid.get_size()
        if (!data.empty())
            memcpy(grid.get(), data.data(), std::min<size_t>(grid.get_size(), data.size()) * sizeof(float));
        Manifold& m_ref = *m_ptr;
        m_ref = volume_polygonize(xform, grid, tau, make_triangles, high_is_inside, dual_connectivity);
    }

    void graph_to_feq(AMGraph3D* g_ptr, Manifold* m_ptr, const std::vector<double>& node_radii, bool symmetrize, bool use_graph_radii) {
        vector<double> node_rs;
        if (use_graph_radii) {
            const size_t N = g_ptr->no_nodes();
            node_rs.resize(N);
            for (auto n : g_ptr->node_ids())
                node_rs[n] = g_ptr->node_color[n][1];
        } else {
            node_rs = node_radii; // copy provided radii
        }

        *m_ptr = graph_to_FEQ(*g_ptr, node_rs, symmetrize);
    }

    void non_rigid_registration(Manifold* m_ptr, Manifold* m_ref_ptr) {
        ::non_rigid_registration(*m_ptr, *m_ref_ptr);
    }

    void rsr_recon(Manifold* m_ptr,
                const std::vector<double>& verts,
                const std::vector<double>& normals,
                int v_num, int n_num,
                bool isEuclidean, int genus, int k, int r, int theta, int n) {
        vector<Vec3d> vertices;
        vector<Vec3d> norm;
        vertices.reserve(v_num);
        norm.reserve(n_num);
        // Expecting packed as [x0,y0,z0,x1,y1,z1,...]
        for (int i = 0; i < v_num; i++) {
            size_t base = static_cast<size_t>(3 * i);
            if (base + 2 < verts.size())
                vertices.emplace_back(verts[base], verts[base + 1], verts[base + 2]);
        }
        for (int i = 0; i < n_num; i++) {
            size_t base = static_cast<size_t>(3 * i);
            if (base + 2 < normals.size())
                norm.emplace_back(normals[base], normals[base + 1], normals[base + 2]);
        }
        reconstruct_single(*m_ptr, vertices, norm,
                        isEuclidean, genus, k, r, theta, n);
    }

    void extrude_faces(Manifold* m_ptr, const std::vector<int>& faces, std::vector<size_t>& fidx_ref) {
        FaceSet fset;
        for (int fid : faces)
            fset.insert(FaceID(fid));
        FaceSet floop = extrude_face_set(*m_ptr, fset);
        for (auto f : floop)
            fidx_ref.push_back(f.index);
    }

    void kill_face_loop(Manifold* m_ptr) {
        kill_face_loop(*m_ptr);
    }

    void kill_degenerate_face_loops(Manifold* m_ptr, double thresh) {
        kill_degenerate_face_loops(*m_ptr, thresh);
    }

    void stable_marriage_registration(Manifold* m_ptr, Manifold* m_ref_ptr) {
        ::stable_marriage_registration(*m_ptr, *m_ref_ptr);
    }

    vector<Manifold*>* connected_components(Manifold* m_ptr) {
        vector<Manifold*>* components = new vector<Manifold*>();
        // Compute connected components and fill the components vector
        for (auto comp: connected_components(*m_ptr)) {
            components->push_back(new Manifold(comp));
        }
        return components;
    }

    size_t mesh_vec_size(const vector<Manifold*>* mv_ptr) {
        return mv_ptr->size();
    }

    Manifold* mesh_vec_get(const vector<Manifold*>* mv_ptr, size_t i) {
        return mv_ptr->at(i);
    }

    void mesh_vec_del(vector<Manifold*>* mv_ptr) {
        delete mv_ptr;
    }

    int count_boundary_curves(Manifold* m_ptr) {
        return count_boundary_curves(*m_ptr);
    }
}
