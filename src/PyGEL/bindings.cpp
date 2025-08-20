#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>

#include "Graph.h"
#include "Manifold.h"
#include "hmesh_functions.h"
#include "graph_functions.h"
#include "Vec3dVector.h"

#include "I3DTree.h"
#include "MeshDistance.h"
#include "Viewer.h"

namespace py = pybind11;
using namespace PyGEL;
using namespace Geometry;
using namespace HMesh;
using namespace CGLA;

PYBIND11_MODULE(PyGEL, m) {
    m.doc() = "Python bindings for GEL (Geometry and Linear algebra library) - Complete";
    
    // Register opaque types for Python compatibility
    py::class_<Geometry::AMGraph3D>(m, "Graph");
    py::class_<HMesh::Manifold>(m, "Manifold");
    py::class_<GLManifoldViewer>(m, "GLManifoldViewer");
    py::class_<std::vector<CGLA::Vec3d>>(m, "Vec3dVector");
    py::class_<std::vector<HMesh::Manifold*>>(m, "MeshVec");
    
    // Constants
    m.attr("InvalidIndex") = PyGEL::InvalidIndex;
    
    // GRAPH FUNCTIONS - Using existing wrapper functions
    m.def("Graph_new", []() -> AMGraph3D* {
        return Graph_new();
    });
    m.def("Graph_delete", [](AMGraph3D* self) {
        Graph_delete(self);
    });
    m.def("Graph_clear", [](AMGraph3D* self) {
        Graph_clear(self);
    });
    m.def("Graph_add_node", [](AMGraph3D* self, const std::vector<double>& pos) -> size_t {
        return Graph_add_node(self, pos);
    });
    m.def("Graph_connect_nodes", [](AMGraph3D* self, size_t n0, size_t n1) -> size_t {
        return Graph_connect_nodes(self, n0, n1);
    });
    m.def("Graph_disconnect_nodes", [](AMGraph3D* self, size_t n0, size_t n1) {
        Graph_disconnect_nodes(self, n0, n1);
    });
    m.def("Graph_remove_node", [](AMGraph3D* self, size_t n) {
        Graph_remove_node(self, n);
    });
    m.def("Graph_nodes", [](AMGraph3D* self) -> std::vector<size_t> {
        return Graph_nodes(self);
    });
    m.def("Graph_neighbors", [](AMGraph3D* self, size_t n, const std::string& mode) -> std::vector<size_t> {
        char mchar = mode.empty() ? 'n' : mode[0];
        return Graph_neighbors(self, n, mchar);
    });
    m.def("Graph_node_in_use", [](AMGraph3D* self, size_t n) -> bool {
        return Graph_node_in_use(self, n);
    });
    m.def("Graph_average_edge_length", [](AMGraph3D* self) -> double {
        return Graph_average_edge_length(self);
    });
    m.def("Graph_copy", [](AMGraph3D* self) -> AMGraph3D* {
        return Graph_copy(self);
    });
    m.def("Graph_merge_nodes", [](AMGraph3D* self, size_t n0, size_t n1, bool avg) {
        Graph_merge_nodes(self, n0, n1, avg);
    });
    m.def("Graph_no_nodes", [](AMGraph3D* self) -> size_t {
        return self->no_nodes();
    });
    m.def("Graph_no_edges", [](AMGraph3D* self) -> size_t {
        return self->no_edges();
    });

    // MANIFOLD FUNCTIONS - Using existing wrapper functions where available
    m.def("Manifold_new", []() -> Manifold* {
        return Manifold_new();
    });
    m.def("Manifold_delete", [](Manifold* self) {
        Manifold_delete(self);
    });
    m.def("Manifold_copy", [](Manifold* self) -> Manifold* {
        return Manifold_copy(self);
    });
    m.def("Manifold_from_triangles", [](const std::vector<double>& vertices, const std::vector<int>& faces) -> Manifold* {
        return Manifold_from_triangles(vertices, faces);
    });
    m.def("Manifold_from_points", [](int N, const std::vector<double>& pts, const std::vector<double>& X_axis, const std::vector<double>& Y_axis) -> Manifold* {
        return Manifold_from_points(N, pts, Vec(X_axis[0], X_axis[1], X_axis[2]), Vec(Y_axis[0], Y_axis[1], Y_axis[2]));
    });
    m.def("Manifold_merge", [](Manifold* self, Manifold* other) {
        Manifold_merge(self, other);
    });
    
    // ...existing code...
    m.def("Manifold_positions", [](Manifold* self) -> py::array_t<Scalar> {
        return Manifold_positions(self);
    });
    
    m.def("Manifold_vertices", [](Manifold* self) -> std::vector<size_t> {
        return Manifold_vertices(self);
    });

    m.def("Manifold_faces", [](Manifold* self) -> std::vector<size_t> {
        return Manifold_faces(self);
    });

    m.def("Manifold_halfedges", [](Manifold* self) -> std::vector<size_t> {
        return Manifold_halfedges(self);
    });

    m.def("Manifold_circulate_vertex", [](Manifold* self, size_t v, const std::string& mode) -> std::vector<size_t> {
        char mchar = mode.empty() ? 'v' : mode[0];
        return Manifold_circulate_vertex(self, v, mchar);
    });

    m.def("Manifold_circulate_face", [](Manifold* self, size_t f, const std::string& mode) -> std::vector<size_t> {
        char mchar = mode.empty() ? 'v' : mode[0];
        return Manifold_circulate_face(self, f, mchar);
    });
    
    m.def("Manifold_add_face", [](Manifold* self, py::buffer pos_buf) -> size_t {
        py::buffer_info info = pos_buf.request();
        
        // Ensure it's a double buffer
        if (info.format != py::format_descriptor<double>::format()) {
            throw std::runtime_error("Incompatible buffer format: expected double array");
        }
        
        // Create vector from buffer
        std::vector<double> pos(static_cast<double*>(info.ptr), 
                               static_cast<double*>(info.ptr) + info.size);
        
        return Manifold_add_face(self, pos);
    });
    
    m.def("Manifold_split_face_by_edge", [](Manifold* self, size_t f, size_t v0, size_t v1) -> size_t {
        return Manifold_split_face_by_edge(self, f, v0, v1);
    });
    
    m.def("Manifold_split_face_by_vertex", [](Manifold* self, size_t f) -> size_t {
        return Manifold_split_face_by_vertex(self, f);
    });
    
    m.def("Manifold_stitch_boundary_edges", [](Manifold* self, size_t h0, size_t h1) -> bool {
        return Manifold_stitch_boundary_edges(self, h0, h1);
    });
    
    m.def("Manifold_merge_faces", [](Manifold* self, size_t f, size_t h) -> bool {
        return Manifold_merge_faces(self, f, h);
    });
    
    // Additional geometry functions
    m.def("vertex_normal", [](Manifold* self, size_t v) -> HMesh::Manifold::Vec {
        return vertex_normal(self, v);
    });
    
    m.def("face_normal", [](Manifold* self, size_t f) -> HMesh::Manifold::Vec {
        return face_normal(self, f);
    });
    
    m.def("principal_curvatures", [](Manifold* self, size_t v) -> std::vector<double> {
        return principal_curvatures(self, v);
    });
    
    m.def("centre", [](Manifold* self, size_t f) -> HMesh::Manifold::Vec {
        return centre(self, f);
    });
    
    m.def("Manifold_no_vertices", [](Manifold* self) -> size_t {
        return self->no_vertices();
    });
    
    m.def("Manifold_no_faces", [](Manifold* self) -> size_t {
        return self->no_faces();
    });
    
    m.def("Manifold_no_halfedges", [](Manifold* self) -> size_t {
        return self->no_halfedges();
    });
    
    m.def("Manifold_allocated_vertices", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_vertices(self);
    });
    
    m.def("Manifold_allocated_faces", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_faces(self);
    });
    
    m.def("Manifold_allocated_halfedges", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_halfedges(self);
    });
    
    m.def("Manifold_clear", [](Manifold* self) {
        self->clear();
    });
    
    m.def("Manifold_cleanup", [](Manifold* self) {
        self->cleanup();
    });
    
    // Additional Manifold operations
    m.def("Manifold_close_hole", [](Manifold* self, size_t h) -> size_t {
        return Manifold_close_hole(self, h);
    });
    
    // Manifold editing functions - using wrapper functions
    m.def("Manifold_remove_face", [](Manifold* self, size_t fid) -> bool {
        return Manifold_remove_face(self, fid);
    });
    
    m.def("Manifold_remove_vertex", [](Manifold* self, size_t vid) -> bool {
        return Manifold_remove_vertex(self, vid);
    });
    
    m.def("Manifold_remove_edge", [](Manifold* self, size_t hid) -> bool {
        return Manifold_remove_edge(self, hid);
    });
    
    m.def("Manifold_vertex_in_use", [](Manifold* self, size_t vid) -> bool {
        return Manifold_vertex_in_use(self, vid);
    });
    
    m.def("Manifold_face_in_use", [](Manifold* self, size_t fid) -> bool {
        return Manifold_face_in_use(self, fid);
    });
    
    m.def("Manifold_halfedge_in_use", [](Manifold* self, size_t hid) -> bool {
        return Manifold_halfedge_in_use(self, hid);
    });
    
    // Basic topology queries - using wrapper functions where available
    m.def("Manifold_is_vertex_at_boundary", [](Manifold* self, size_t vid) -> bool {
        return is_vertex_at_boundary(self, vid);
    });
    
    m.def("Manifold_is_halfedge_at_boundary", [](Manifold* self, size_t hid) -> bool {
        return is_halfedge_at_boundary(self, hid);
    });
    
    m.def("Manifold_connected", [](Manifold* self, size_t v0, size_t v1) -> bool {
        return connected(self, v0, v1);
    });
    
    // Geometry and measurement functions
    m.def("length", [](Manifold* self, size_t h) -> double {
        return length(self, h);
    });
    
    m.def("valency", [](Manifold* self, size_t v) -> size_t {
        return valency(self, v);
    });
    
    m.def("boundary_edge", [](Manifold* self, size_t v, size_t h) -> bool {
        return boundary_edge(self, v, h);
    });
    
    m.def("no_edges", [](Manifold* self, size_t f) -> size_t {
        return no_edges(self, f);
    });
    
    m.def("area", [](Manifold* self, size_t f) -> double {
        return area(self, f);
    });
    
    m.def("one_ring_area", [](Manifold* self, size_t v) -> double {
        return one_ring_area(self, v);
    });
    
    m.def("mixed_area", [](Manifold* self, size_t v) -> double {
        return mixed_area(self, v);
    });
    
    m.def("gaussian_curvature", [](Manifold* self, size_t v) -> double {
        return gaussian_curvature(self, v);
    });
    
    m.def("mean_curvature", [](Manifold* self, size_t v) -> double {
        return mean_curvature(self, v);
    });
    
    m.def("perimeter", [](Manifold* self, size_t f) -> double {
        return perimeter(self, f);
    });
    
    m.def("total_area", [](Manifold* self) -> double {
        return total_area(self);
    });
    
    m.def("volume", [](Manifold* self) -> double {
        return volume(self);
    });
    
    // Walker functions - using wrapper functions
    m.def("Walker_next_halfedge", [](Manifold* self, size_t h) -> size_t {
        return Walker_next_halfedge(self, h);
    });
    
    m.def("Walker_prev_halfedge", [](Manifold* self, size_t h) -> size_t {
        return Walker_prev_halfedge(self, h);
    });
    
    m.def("Walker_opposite_halfedge", [](Manifold* self, size_t h) -> size_t {
        return Walker_opposite_halfedge(self, h);
    });
    
    m.def("Walker_incident_face", [](Manifold* self, size_t h) -> size_t {
        return Walker_incident_face(self, h);
    });
    
    m.def("Walker_incident_vertex", [](Manifold* self, size_t h) -> size_t {
        return Walker_incident_vertex(self, h);
    });
    
    // Mesh modification functions - using wrapper functions
    m.def("Manifold_flip_edge", [](Manifold* self, size_t h) -> bool {
        return Manifold_flip_edge(self, h);
    });
    
    m.def("Manifold_collapse_edge", [](Manifold* self, size_t h, bool avg_vertices) -> bool {
        return Manifold_collapse_edge(self, h, avg_vertices);
    });
    
    m.def("Manifold_split_edge", [](Manifold* self, size_t h) -> size_t {
        return Manifold_split_edge(self, h);
    });
    
    // HMESH FUNCTIONS - Basic mesh validation and processing
    m.def("valid", [](Manifold* self) -> bool {
        return valid(self);
    });
    
    m.def("closed", [](Manifold* self) -> bool {
        return closed(self);
    });
    
    // Additional HMesh processing functions
    m.def("stitch_mesh", [](Manifold* self, double rad) -> int {
        return stitch_mesh(self, rad);
    });
    
    m.def("count_boundary_curves", [](Manifold* self) -> int {
        return count_boundary_curves(self);
    });
    
    // Mesh processing functions
    m.def("remove_caps", [](Manifold* self, float thresh) {
        remove_caps(self, thresh);
    });
    
    m.def("remove_needles", [](Manifold* self, float thresh, bool avgPos) {
        remove_needles(self, thresh, avgPos);
    });
    
    m.def("close_holes", [](Manifold* self, int max_size) {
        close_holes(self, max_size);
    });
    
    m.def("flip_orientation", [](Manifold* self) {
        flip_orientation(self);
    });
    
    m.def("merge_coincident_boundary_vertices", [](Manifold* self, double rad) {
        merge_coincident_boundary_vertices(self, rad);
    });
    
    // Mesh optimization functions
    m.def("minimize_curvature", [](Manifold* self, bool anneal) {
        minimize_curvature(self, anneal);
    });
    
    m.def("minimize_dihedral_angle", [](Manifold* self, int max_iter, bool anneal, bool alpha, double gamma) {
        minimize_dihedral_angle(self, max_iter, anneal, alpha, gamma);
    });
    
    m.def("maximize_min_angle", [](Manifold* self, float thresh, bool anneal) {
        maximize_min_angle(self, thresh, anneal);
    });
    
    m.def("optimize_valency", [](Manifold* self, bool anneal) {
        optimize_valency(self, anneal);
    });
    
    m.def("randomize_mesh", [](Manifold* self, int max_iter) {
        randomize_mesh(self, max_iter);
    });
    
    // Mesh simplification and refinement
    m.def("quadric_simplify", [](Manifold* self, double keep_fraction, double singular_thresh, double error_thresh) {
        quadric_simplify(self, keep_fraction, singular_thresh, error_thresh);
    });
    
    m.def("average_edge_length", [](Manifold* self) -> float {
        return average_edge_length(self);
    });
    
    m.def("median_edge_length", [](Manifold* self) -> float {
        return median_edge_length(self);
    });
    
    m.def("refine_edges", [](Manifold* self, float t) -> int {
        return refine_edges(self, t);
    });
    
    // Subdivision functions
    m.def("cc_split", [](Manifold* self) {
        cc_split(self);
    });
    
    m.def("loop_split", [](Manifold* self) {
        loop_split(self);
    });
    
    m.def("root3_subdivide", [](Manifold* self) {
        root3_subdivide(self);
    });
    
    m.def("rootCC_subdivide", [](Manifold* self) {
        rootCC_subdivide(self);
    });
    
    m.def("butterfly_subdivide", [](Manifold* self) {
        butterfly_subdivide(self);
    });
    
    // Smoothing functions
    m.def("cc_smooth", [](Manifold* self) {
        cc_smooth(self);
    });
    
    m.def("volume_preserving_cc_smooth", [](Manifold* self, int iter) {
        volume_preserving_cc_smooth(self, iter);
    });
    
    m.def("loop_smooth", [](Manifold* self) {
        loop_smooth(self);
    });
    
    m.def("taubin_smooth", [](Manifold* self, int iter) {
        taubin_smooth(self, iter);
    });
    
    m.def("laplacian_smooth", [](Manifold* self, float weight, int iter) {
        laplacian_smooth(self, weight, iter);
    });
    
    m.def("anisotropic_smooth", [](Manifold* self, float sharpness, int iter) {
        anisotropic_smooth(self, sharpness, iter);
    });
    
    // Triangulation functions
    m.def("shortest_edge_triangulate", [](Manifold* self) {
        shortest_edge_triangulate(self);
    });
    
    m.def("ear_clip_triangulate", [](Manifold* self) {
        ear_clip_triangulate(self);
    });
    
    // GRAPH FUNCTIONS - Graph processing
    m.def("graph_from_mesh", [](Manifold* m, AMGraph3D* g) {
        graph_from_mesh(m, g);
    });
    
    m.def("graph_load", [](AMGraph3D* g, const std::string& filename) -> bool {
        return graph_load(g, filename);
    });
    
    m.def("graph_save", [](AMGraph3D* g, const std::string& filename) -> bool {
        return graph_save(g, filename);
    });
    
    m.def("graph_to_mesh_cyl", [](AMGraph3D* g, Manifold* m, float fudge) {
        graph_to_mesh_cyl(g, m, fudge);
    });
    
    m.def("graph_smooth", [](AMGraph3D* g, int iter, float alpha) {
        graph_smooth(g, iter, alpha);
    });
    
    m.def("graph_edge_contract", [](AMGraph3D* g, double dist_thresh) -> int {
        return graph_edge_contract(g, dist_thresh);
    });
    
    m.def("graph_prune", [](AMGraph3D* g) {
        graph_prune(g);
    });
    
    m.def("graph_saturate", [](AMGraph3D* g, int hops, double dist_frac, double rad) {
        graph_saturate(g, hops, dist_frac, rad);
    });
    
    // VECTOR FUNCTIONS - Vec3dVector utilities only
    
    m.def("Vec3dVector_new", [](size_t s) -> Vec3dVector* {
        return Vec3dVector_new(s);
    });
    
    m.def("Vec3dVector_size", [](Vec3dVector* self) -> size_t {
        return Vec3dVector_size(self);
    });
    
    m.def("Vec3dVector_delete", [](Vec3dVector* self) {
        Vec3dVector_delete(self);
    });
    
    // Additional Manifold functions with wrapper functions
    m.def("Manifold_no_allocated_vertices", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_vertices(self);
    });
    
    m.def("Manifold_no_allocated_faces", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_faces(self);
    });
    
    m.def("Manifold_no_allocated_halfedges", [](Manifold* self) -> size_t {
        return Manifold_no_allocated_halfedges(self);
    });
    
    // MISSING HMESH FUNCTIONS
    
    // File I/O functions
    m.def("load", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return load(filename, m_ptr);
    });
    
    m.def("obj_load", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return obj_load(filename, m_ptr);
    });
    
    m.def("off_load", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return off_load(filename, m_ptr);
    });
    
    m.def("ply_load", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return ply_load(filename, m_ptr);
    });
    
    m.def("x3d_load", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return x3d_load(filename, m_ptr);
    });
    
    m.def("obj_save", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return obj_save(filename, m_ptr);
    });
    
    m.def("off_save", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return off_save(filename, m_ptr);
    });
    
    m.def("x3d_save", [](const std::string& filename, Manifold* m_ptr) -> bool {
        return x3d_save(filename, m_ptr);
    });
    
    // Bounding box and sphere functions
    m.def("bbox", [](Manifold* m_ptr) -> std::pair<std::vector<double>, std::vector<double>> {
        return bbox(m_ptr);
    });
    
    m.def("bsphere", [](Manifold* m_ptr) -> std::pair<std::vector<double>, double> {
        return bsphere(m_ptr);
    });
    
    // Additional smoothing function
    m.def("regularize_quads", [](Manifold* self, float weight, float shrink, int iter) {
        regularize_quads(self, weight, shrink, iter);
    });
    
    // Volumetric isocontour function
    m.def("volumetric_isocontour", [](Manifold* m_ptr, int x_dim, int y_dim, int z_dim, 
                                      const std::vector<float>& data, const std::vector<double>& pmin, 
                                      const std::vector<double>& pmax, float tau, bool make_triangles, 
                                      bool high_is_inside, bool dual_connectivity) {
        volumetric_isocontour(m_ptr, x_dim, y_dim, z_dim, data, pmin, pmax, tau, make_triangles, high_is_inside, dual_connectivity);
    });
    
    // Graph to mesh conversion
    m.def("graph_to_feq", [](AMGraph3D* g_ptr, Manifold* m_ptr, const std::vector<double>& node_radii, 
                              bool symmetrize, bool use_graph_radii) {
        graph_to_feq(g_ptr, m_ptr, node_radii, symmetrize, use_graph_radii);
    });
    
    // Mesh registration and reconstruction
    m.def("non_rigid_registration", [](Manifold* m_ptr, Manifold* m_ref_ptr) {
        non_rigid_registration(m_ptr, m_ref_ptr);
    });
    
    m.def("rsr_recon", [](Manifold* m_ptr, const std::vector<double>& verts, const std::vector<double>& normals, 
                          int v_num, int n_num, bool isEuclidean, int genus, int k, int r, int theta, int n) {
        rsr_recon(m_ptr, verts, normals, v_num, n_num, isEuclidean, genus, k, r, theta, n);
    });
    
    // Face extrusion and loop operations
    m.def("extrude_faces", [](Manifold* m_ptr, const std::vector<int>& faces, std::vector<size_t>& fidx_ptr) {
        extrude_faces(m_ptr, faces, fidx_ptr);
    });
    
    m.def("kill_face_loop", [](Manifold* m_ptr) {
        kill_face_loop(m_ptr);
    });
    
    m.def("kill_degenerate_face_loops", [](Manifold* m_ptr, double thresh) {
        kill_degenerate_face_loops(m_ptr, thresh);
    });
    
    m.def("stable_marriage_registration", [](Manifold* m_ptr, Manifold* m_ref_ptr) {
        stable_marriage_registration(m_ptr, m_ref_ptr);
    });
    
    // Connected components
    m.def("connected_components", [](Manifold* m_ptr) -> std::vector<Manifold*> {
        return connected_components(*m_ptr);
    });
    
    m.def("mesh_vec_size", [](std::vector<Manifold*>* mv_ptr) -> size_t {
        return mesh_vec_size(*mv_ptr);
    });
    
    m.def("mesh_vec_get", [](std::vector<Manifold*>* mv_ptr, size_t i) -> Manifold* {
        return mesh_vec_get(*mv_ptr, i);
    });
    
    m.def("mesh_vec_del", [](std::vector<Manifold*>* mv_ptr) {
        mesh_vec_del(mv_ptr);
    });
    

    // I3DTree bindings
    py::class_<I3DTree>(m, "I3DTree");
    m.def("I3DTree_new", []() -> I3DTree* {
        return I3DTree_new();
    });
    m.def("I3DTree_delete", [](I3DTree* self) {
        I3DTree_delete(self);
    });
    m.def("I3DTree_insert", [](I3DTree* tree, double x, double y, double z, size_t v) {
        I3DTree_insert(tree, x, y, z, v);
    });
    m.def("I3DTree_build", [](I3DTree* tree) {
        I3DTree_build(tree);
    });
    m.def("I3DTree_closest_point", [](I3DTree* tree, double x, double y, double z, double r) -> std::pair<std::vector<double>, size_t> {
        return I3DTree_closest_point(tree, x, y, z, r);
    });
    m.def("I3DTree_in_sphere", [](I3DTree* tree, double x, double y, double z, double r) -> std::pair<std::vector<CGLA::Vec3d>, std::vector<size_t>> {
        return I3DTree_in_sphere(tree, x, y, z, r);
    });
    m.def("I3DTree_m_closest_points", [](I3DTree* tree, double x, double y, double z, double r, int m) -> std::pair<std::vector<CGLA::Vec3d>, std::vector<size_t>> {
        return I3DTree_m_closest_points(tree, x, y, z, r, m);
    });
    
    // Graph to mesh functions
    m.def("graph_to_mesh_iso", [](AMGraph3D* g_ptr, Manifold* m_ptr, float fudge, size_t grid_res) {
        graph_to_mesh_iso(g_ptr, m_ptr, fudge, grid_res);
    });
    
    // Skeletonization functions
    m.def("graph_LS_skeleton", [](AMGraph3D* g_ptr, AMGraph3D* skel_ptr, bool sampling) -> std::vector<size_t> {
        return graph_LS_skeleton(g_ptr, skel_ptr, sampling);
    });
    m.def("graph_MSLS_skeleton", [](AMGraph3D* g_ptr, AMGraph3D* skel_ptr, int grow_thresh) -> std::vector<size_t> {
        return graph_MSLS_skeleton(g_ptr, skel_ptr, grow_thresh);
    });
    m.def("graph_front_skeleton", [](AMGraph3D* g_ptr, AMGraph3D* skel_ptr, int N_col, const std::vector<double>& colors, int intervals) -> std::vector<size_t> {
        return graph_front_skeleton(g_ptr, skel_ptr, N_col, colors, intervals);
    });
    m.def("graph_combined_skeleton", [](AMGraph3D* g_ptr, AMGraph3D* skel_ptr, int N_col, const std::vector<double>& colors, int intervals) -> std::vector<size_t> {
        return graph_combined_skeleton(g_ptr, skel_ptr, N_col, colors, intervals);
    });
    
    m.def("graph_minimum_spanning_tree", [](AMGraph3D* g_ptr, AMGraph3D* mst_ptr, int root) {
        graph_minimum_spanning_tree(g_ptr, mst_ptr, root);
    });
    
    m.def("graph_close_chordless_cycles", [](AMGraph3D* g_ptr, int root, int hops, double rad) {
        graph_close_chordless_cycles(g_ptr, root, hops, rad);
    });
    
    // MISSING MESHDISTANCE FUNCTIONS
    
    // MeshDistance class functions - using opaque pointers
    m.def("MeshDistance_new", [](Manifold* m) -> MeshDistance* {
        return MeshDistance_new(m);
    });
    
    m.def("MeshDistance_delete", [](MeshDistance* self) {
        MeshDistance_delete(self);
    });
    
    m.def("MeshDistance_signed_distance", [](MeshDistance* self, const std::vector<float>& p, float upper) -> std::vector<float> {
        return MeshDistance_signed_distance(self, p, upper);
    });
    
    m.def("MeshDistance_ray_inside_test", [](MeshDistance* self, const std::vector<float>& p, int no_rays) -> std::vector<int> {
        return MeshDistance_ray_inside_test(self, p, no_rays);
    });
    
    m.def("MeshDistance_ray_intersect", [](MeshDistance* self, const std::vector<float>& _p, const std::vector<float>& _d) -> std::pair<bool, float> {
        CGLA::Vec3f p(_p[0], _p[1], _p[2]);
        CGLA::Vec3f d(_d[0], _d[1], _d[2]);
        float t;
        bool result = MeshDistance_ray_intersect(self, p, d, &t);
        return {result, t};
    });
    
    // Test functions
    m.def("test_complete", []() -> std::string {
        return "PyGEL complete version working!";
    });
    
    m.def("test_graph", []() -> std::string {
        // Test graph creation using wrapper functions
        auto g = Graph_new();
        auto n1 = Graph_add_node(g, {1, 2, 3});
        auto n2 = Graph_add_node(g, {4, 5, 6});
        Graph_connect_nodes(g, n1, n2);
        auto result = "Graph test: " + std::to_string(g->no_nodes()) + " nodes, " + std::to_string(g->no_edges()) + " edges";
        Graph_delete(g);
        return result;
    });
    
    m.def("test_manifold", []() -> std::string {
        // Test manifold creation using wrapper functions
        auto m = Manifold_new();
        auto result = "Manifold test: " + std::to_string(m->no_vertices()) + " vertices, " + std::to_string(m->no_faces()) + " faces";
        Manifold_delete(m);
        return result;
    });

    // VIEWER (GLManifoldViewer) bindings (only if built with GLGraphics)
    m.def("GLManifoldViewer_new", []() -> GLManifoldViewer* { return GLManifoldViewer_new(); });
    m.def("GLManifoldViewer_delete", [](GLManifoldViewer* v) { GLManifoldViewer_delete(v); });
    m.def("GLManifoldViewer_clone_controller", [](GLManifoldViewer* self, GLManifoldViewer* other) {
        GLManifoldViewer_clone_controller(self, other);
    });
    m.def("GLManifoldViewer_display", [](GLManifoldViewer* self,
                                          Manifold* m_ptr,
                                          AMGraph3D* g_ptr,
                                          const std::string& mode,
                                          bool smooth_shading,
                                          const std::vector<float>& bg_color,
                                          std::vector<double>& attrib_vec,
                                          bool reset_view,
                                          bool once) {
        char c = mode.empty() ? 'n' : mode[0];
        Vec3f bg(bg_color[0], bg_color[1], bg_color[2]);
        GLManifoldViewer_display(self, m_ptr, g_ptr, c, smooth_shading, bg, attrib_vec, reset_view, once);
    }, py::arg("viewer"), py::arg("manifold")=nullptr, py::arg("graph")=nullptr, py::arg("mode")="w",
       py::arg("smooth_shading")=true, py::arg("bg_color")=std::vector<float>{0.1f,0.1f,0.1f},
       py::arg("attrib_vec")=std::vector<double>{}, py::arg("reset_view")=false, py::arg("once")=false);

    // m.def("GLManifoldViewer_get_annotation_points", [](GLManifoldViewer_ptr self) -> std::vector<double> {
    //     return GLManifoldViewer_get_annotation_points(self);
    // });

    // m.def("GLManifoldViewer_set_annotation_points", [](GLManifoldViewer_ptr self, const std::vector<double>& data) {
    //     GLManifoldViewer_set_annotation_points(self, data);
    // });
}
