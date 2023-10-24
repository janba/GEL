#include <array>
#include <cmath>
#include <algorithm>

#include <GEL/CGLA/CGLA.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/SphereDelaunay.h>
#include <GEL/HMesh/comb_quad.h>

using namespace Geometry;
using namespace CGLA;
using namespace HMesh;
using namespace std;
using NodeID = AMGraph::NodeID;
using NodeSet = AMGraph::NodeSet;

// Initialize global arrays

map<NodeID, int> val2deg;
map<pair<NodeID,NodeID>, int> branchdeg;
map<pair<NodeID,NodeID>, HMesh::FaceID> branchface;
map<pair<NodeID,NodeID>, HMesh::FaceID> branch_best_face;
map<pair<NodeID,NodeID>, HMesh::VertexID> branch_best_vertex;
map<pair<NodeID,NodeID>, CGLA::Vec3d> branch2vert;
map<FaceID, VertexID> one_ring_face_vertex;

void clear_global_arrays() {
    val2deg.clear();
    branchdeg.clear();
    branchface.clear();
    branch_best_face.clear();
    branch_best_vertex.clear();
    branch2vert.clear();
    one_ring_face_vertex.clear();
    return;
}

//Graph util functions

vector<NodeID> next_neighbours(const Geometry::AMGraph3D& g, NodeID prev, NodeID curr) {

    vector<NodeID> neighbour_list;
    auto N = g.neighbors(curr);
    for (auto next: N) {
        if(next != prev)
            neighbour_list.push_back(next);
    }
    return neighbour_list;

}


//Mesh util functions

void id_preserving_cc(HMesh::Manifold& m_in) {

    vector<FaceID> base_faces;
    int Invalid = -1;
    HalfEdgeAttributeVector<int> htouched(m_in.allocated_halfedges(), Invalid);
    int Valid = 1;
    map<FaceID,VertexID> face2centerv;
    vector<HalfEdgeID> base_edges;

    vector<HalfEdgeID> new_edges;


    for(auto f: m_in.faces())
        base_faces.push_back(f);

    for(auto h: m_in.halfedges())
        base_edges.push_back(h);

    for(auto f: base_faces)
        if(m_in.in_use(f)) {
            VertexID center_v = m_in.split_face_by_vertex(f);
            for(Walker w = m_in.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                new_edges.push_back(w.halfedge());
                face2centerv.insert(std::make_pair(w.face(), center_v));
            }
        }

    FaceAttributeVector<int> ftouched(m_in.allocated_faces(), Invalid);

    for (auto h: base_edges) {
        FaceID f1 = m_in.walker(h).face();
        FaceID f2 = m_in.walker(h).opp().face();
        if(ftouched[f1]==Invalid && ftouched[f2]==Invalid) {
            VertexID opp_v = m_in.split_edge(h);
            ftouched[f1] = Valid;
            ftouched[f2] = Valid;
            m_in.split_face_by_edge(f1, face2centerv.find(f1)->second, opp_v);
            m_in.split_face_by_edge(f2, face2centerv.find(f2)->second, opp_v);
        }

    }
    for (auto h_dissolve: new_edges)
        if(m_in.in_use(h_dissolve))
            m_in.merge_faces(m_in.walker(h_dissolve).face(), h_dissolve);
    return;

}

void quad_mesh_leaves(HMesh::Manifold& m, VertexAttributeVector<NodeID>& vertex2node) {

    vector<FaceID> base_faces;
    vector<HalfEdgeID> new_edges;

    for(auto f: m.faces())
        if(no_edges(m, f) != 4) {
            HalfEdgeID ref_h;

            VertexID ref_v  = one_ring_face_vertex[f];

            if(ref_v == InvalidVertexID)
                continue;

            VertexID center_v = m.split_face_by_vertex(f);
            vertex2node[center_v] = vertex2node[m.walker(center_v).vertex()];
            int counter = 0;
            int dissolve_flag = 0;
            for(Walker w = m.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                if(m.walker(w.halfedge()).vertex() == ref_v || m.walker(w.halfedge()).opp().vertex() == ref_v) {
                    dissolve_flag = counter%2;
                    ref_h = w.halfedge();
                }
                counter++;
            }
            counter = 0;
            for(Walker w = m.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                if(counter%2 != dissolve_flag)
                    new_edges.push_back(w.halfedge());
                counter++;
            }

        }


    for (auto h_dissolve: new_edges)
        if(m.in_use(h_dissolve))
            m.merge_faces(m.walker(h_dissolve).face(), h_dissolve);


    return;


}



//Graph - Mesh relationship Functions

VertexID branch2vertex (const HMesh::Manifold &m_out, const Geometry::AMGraph3D& g,
                        NodeID n, NodeID nn, const Util::AttribVec<NodeID, FaceSet>& node2fs) {

    Vec3d vert_pos = branch2vert.find(std::make_pair(n,nn))->second;

    for (auto v: m_out.vertices())
        if(sqr_length(m_out.pos(v) - vert_pos) == 0)
            return v;

    return InvalidVertexID;

}

void init_branch_degree(const HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                        const Util::AttribVec<NodeID, FaceSet>& node2fs) {


    for (auto n:g.node_ids()) {
        auto N = g.neighbors(n);

        //for all branch nodes

        if(N.size() > 2) {


            // for each outgoing arc

            for (auto nn: N) {

                int src_branch_degree = valency(m, branch2vertex(m, g, n,nn, node2fs));
                vector<NodeID> branch_path;
                NodeID curr_node = nn;
                NodeID prev_node = n;

                int leaf_flag = 0;

                //traverse val 2 nodes to next branch node


                while(true) {
                    auto curr_nbs = next_neighbours(g, prev_node, curr_node);
                    if(curr_nbs.size() > 1)
                        break;
                    else if(curr_nbs.size() == 0) {
                        branch_path.push_back(curr_node); leaf_flag = 1; break;
                    }
                    else {
                        branch_path.push_back(curr_node);
                        prev_node = curr_node;
                        curr_node = curr_nbs[0];
                    }
                }

                //pick lower degree

                int dest_branch_degree;
                if(leaf_flag == 1) {
                    dest_branch_degree = src_branch_degree;
                }
                else
                    dest_branch_degree = valency(m, branch2vertex(m, g, curr_node, prev_node, node2fs));


                int path_degree = 0;
                int jn_degree = 0;

                if(dest_branch_degree < src_branch_degree) {
                    path_degree = (dest_branch_degree)*2;
                    jn_degree = dest_branch_degree - 1;
                    cout << "dst < src " << endl;
                    cout << "dst " << dest_branch_degree << endl;
                    cout << "src " << src_branch_degree << endl;
                }

                else if (dest_branch_degree == src_branch_degree) {
                    path_degree = dest_branch_degree*2;
                    jn_degree = dest_branch_degree;
                }

                else {
                    jn_degree = src_branch_degree - 1;
                    path_degree = (src_branch_degree)*2;
                    cout << "dst > src " << endl;
                    cout << "dst " << dest_branch_degree << endl;
                    cout << "src " << src_branch_degree << endl;
                }

                auto key = std::make_pair(n,nn);
                branchdeg.insert(std::make_pair(key,jn_degree));
                for (auto val2node : branch_path)
                    val2deg.insert(std::make_pair(val2node,path_degree));

            }
        }
    }

    // for junction-less graphs

    bool has_junction = false;

    for (auto n: g.node_ids()) {
        if(g.valence(n) > 2)
            has_junction = true;
    }

    if(!has_junction)
        for (auto n : g.node_ids())
            if(g.valence(n) <= 2)
                if(val2deg.find(n) == val2deg.end())
                    val2deg.insert(std::make_pair(n,4));

}

FaceID branch2face (const HMesh::Manifold &m_out,
                    const Geometry::AMGraph3D& g, NodeID n, NodeID nn,
                    Util::AttribVec<NodeID, FaceSet>& node2fs) {

    VertexID v = branch2vertex(m_out, g, n, nn, node2fs);
    vector<FaceID> face_set;

    double d_max = FLT_MAX;
    FaceID f_max = InvalidFaceID;
    Vec3d pn = g.pos[n];
    Vec3d pnn = g.pos[nn];
    Vec3d v_n_nn = pnn - pn;


    for(Walker w = m_out.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
        face_set.push_back(w.face());
    }

    for(auto f: face_set) {
        double d = dot(v_n_nn, normal(m_out, f));

        Vec3d face_normal = normal(m_out, f);
        Vec3d face_center = centre(m_out,f);
        float face_plane_d = dot(face_normal, face_center);

        float intersection_x = (face_plane_d - dot(face_normal, pn)) / dot(face_normal, v_n_nn);

        Vec3d intersection_pt = pn + intersection_x*v_n_nn;

        d = sqr_length(face_center - intersection_pt);

        if(d < d_max) {
            f_max = f;
            d_max = d;
        }
    }
    if(g.neighbors(n).size()>2)
        node2fs[n].erase(f_max);
    return f_max;


}

double face_dist(const HMesh::Manifold &m_out, const Geometry::AMGraph3D& g, NodeID n, NodeID nn, FaceID f) {
    Vec3d pn = g.pos[n];
    Vec3d pnn = g.pos[nn];
    Vec3d v_n_nn = pnn - pn;


    if(!m_out.in_use(f))
        return 1;

    Vec3d face_normal = normal(m_out, f);
    Vec3d face_center = centre(m_out,f);
    float face_plane_d = dot(face_normal, face_center);

    float intersection_x = 0;
    if(dot(face_normal, v_n_nn) != 0)
        intersection_x = (face_plane_d - dot(face_normal, pn)) / dot(face_normal, v_n_nn);

    Vec3d intersection_pt = pn + intersection_x*v_n_nn;

    float d = sqr_length(face_center - intersection_pt);

    return d;

}

void init_branch_face_pairs(const HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                            Util::AttribVec<NodeID, FaceSet>& node2fs) {

    for (auto n:g.node_ids()) {

        auto N = g.neighbors(n);

        //for all branch nodes

        if(N.size() > 2) {

            // for each outgoing arc

            for (auto nn: N) {

                auto key = std::make_pair(n,nn);

                FaceID f = branch2face(m, g, n, nn, node2fs);
                branch_best_face.insert(std::make_pair(key,f));

                VertexID v = branch2vertex(m, g, n, nn, node2fs);
                branch_best_vertex.insert(std::make_pair(key,v));
            }
        }

    }
}

//Functions for constructing / editing mesh elements from skeletal nodes

vector<Vec3d> get_face_points(int n) {

    vector<Vec3d> face_vertices;
    double h = 0.5;
    double angle = 0;

    for(int i =0; i<n; i++) {

        Vec3d face_vertex = Vec3d(0, h*cos(angle), h*sin(angle));
        face_vertices.push_back(face_vertex);
        angle+=2*22.0/(7.0*n);
    }
    return face_vertices;
}

vector<FaceID> create_face_pair(Manifold& m, const Vec3d& pos, const Mat3x3d& _R, int axis, int num_sides) {
    if(num_sides == 0) {
        vector<FaceID>fvec;
        return fvec;
    }

    vector<Vec3d> face_points = get_face_points(num_sides);
    Mat3x3d R = _R;
    double det = determinant(R);
    if(abs(det))
        if(det<0) {
            Mat3x3d M = identity_Mat3x3d();
            M[2][2] = -1;
            R = R * M;
        }

    vector<FaceID> fvec;
    vector<Vec3d> front_pts;
    for(int i = 0; i < num_sides; i++)
    {
        Vec3d _p = face_points[i];
        Vec3d p(0);
        p[(0+axis)%3] += _p[0];
        p[(1+axis)%3] += _p[1];
        p[(2+axis)%3] += _p[2];
        front_pts.push_back(R*p+pos);
    }
    fvec.push_back(m.add_face(front_pts));

    vector<Vec3d> back_pts;
    for (int i = 0; i < num_sides; i++) {
        int curr_index = 1 - i;
        if(curr_index < 0 )
            curr_index = num_sides + curr_index;
        Vec3d _p = face_points[curr_index];
        Vec3d p(0);
        p[(0+axis)%3] += _p[0];
        p[(1+axis)%3] += _p[1];
        p[(2+axis)%3] += _p[2];
        back_pts.push_back(R*p+pos);
    }


    fvec.push_back(m.add_face(back_pts));

    for (auto f: fvec) 
        one_ring_face_vertex[f] = m.walker(fvec[0]).vertex();

    return fvec;
}

void val2nodes_to_face_pairs(const Geometry::AMGraph3D &g, HMesh::Manifold &mani,
                             Util::AttribVec<NodeID, FaceSet> &n2fs,
                             VertexAttributeVector<NodeID> &vertex2node,
                             const vector<double> &r)
{
    Util::AttribVec<NodeID, int> touched(g.no_nodes(), 0);
    Util::AttribVec<NodeID, Mat3x3d> warp_frame(g.no_nodes(), identity_Mat3x3d());

    queue<NodeID> Q;

    for (auto middle_node : g.node_ids())
        if (!touched[middle_node])
        {
            Q.push(middle_node);
            while (!Q.empty())
            {
                NodeID n = Q.front();
                Q.pop();
                for (auto m : g.neighbors(n))
                    if (!touched[m])
                    {
                        Q.push(m);
                        touched[m] = 1;
                        Vec3d vect = normalize(g.pos[m] - g.pos[n]);
                        if (g.neighbors(m).size() <= 2)
                        {
                            auto node_list = next_neighbours(g, n, m);
                            if (!node_list.empty())
                            {
                                Vec3d nb_v(0);
                                for (auto m_nb : node_list)
                                    nb_v = g.pos[m_nb];
                                nb_v /= node_list.size();
                                vect += normalize(nb_v - g.pos[m]);
                            }
                        }
                        vect = normalize(vect);
                        Vec3d warp_v = warp_frame[n] * vect;

                        double max_sgn = sign(warp_v[0]);
                        double max_val = abs(warp_v[0]);
                        int max_idx = 0;
                        for (int i = 1; i < 3; ++i)
                        {
                            if (abs(warp_v[i]) > max_val)
                            {
                                max_sgn = sign(warp_v[i]);
                                max_val = abs(warp_v[i]);
                                max_idx = i;
                            }
                        }
                        auto v_target = max_sgn * vect;
                        Quatd q;
                        q.make_rot((warp_frame[n])[max_idx], v_target);

                        warp_frame[m] = transpose(q.get_Mat3x3d() * transpose(warp_frame[n]));

                        if (g.neighbors(m).size() <= 2)
                        {
                            Vec3d s(r[m]);
                            Mat3x3d S = scaling_Mat3x3d(s);
                            auto face_list = create_face_pair(mani, g.pos[m], transpose(warp_frame[m]) * S, max_idx, 8);
                            stitch_mesh(mani, 1e-10);
                            for (auto f : face_list)
                            {
                                n2fs[m].insert(f);
                                for (auto v : mani.incident_vertices(f))
                                    vertex2node[v] = m;
                            }
                        }
                    }
            }
        }
}

int add_ghosts(const vector<Vec3i> &tris, vector<Vec3d> &pts, double thresh)
{

    /* This function creates extra points to add to the BNP vertices for a branch node.
     These extra points are called ghost points because they do not correspond to an outgoing
     edge.
     */

    vector<Vec3d> ghost_pts;

    for (auto t : tris)
    {
        Vec3d v1 = pts[t[1]] - pts[t[0]];
        Vec3d v2 = pts[t[2]] - pts[t[0]];
        Vec3d n = cross(v1, v2);
        ghost_pts.push_back(n);
    }

    /* Next, we cluster the ghost points. This is because in flatish
     configurations we could have several quite similar ghost points. */
    vector<int> cluster_id(ghost_pts.size(), -1);
    int max_id = 0;
    for (int i = 0; i < ghost_pts.size(); ++i)
    {
        if (cluster_id[i] == -1)
        {
            cluster_id[i] = max_id++;
        }
        for (int j = i + 1; j < ghost_pts.size(); ++j)
        {
            if (cluster_id[j] == -1)
            {
                if (dot(normalize(ghost_pts[i]), normalize(ghost_pts[j])) > thresh)
                    cluster_id[j] = cluster_id[i];
            }
        }
    }

    vector<Vec3d> ghost_pts_new(max_id, Vec3d(0));
    for (int i = 0; i < ghost_pts.size(); ++i)
    {
        ghost_pts_new[cluster_id[i]] += ghost_pts[i];
    }

    /* Finally, we cull ghost points too close to an existing non-ghost point.
     The threshold of 0.4 allows ghost points to be a little closer to other
     points than each other. */
    ghost_pts.resize(0);
    for (auto &p : ghost_pts_new)
    {
        p.normalize();
        vector<double> dots;
        for (const auto &p_orig : pts)
            dots.push_back(dot(p, p_orig));
        if (*max_element(begin(dots), end(dots)) < thresh)
            ghost_pts.push_back(p);
    }

    for (auto g : ghost_pts)
        pts.push_back(g);

    return ghost_pts.size();
}


vector<Vec3i> points_to_octahedron(vector<Vec3d> &pts, int s_i, int s_j)
{
    assert(pts.size() == 5);
    Vec3i tri(-1);

    for (int i = 0; i < 5; i++)
    {
        if (i != s_i && i != s_j)
        {
            if (tri[0] == -1)
                tri[0] = i;
            else if (tri[1] == -1)
                tri[1] = i;
            else if (tri[2] == -1)
                tri[2] = i;
        }
    }

    Vec3d n = normalize(cross(pts[tri[1]] - pts[tri[0]], pts[tri[2]] - pts[tri[0]]));

    if (dot(n, pts[s_i]) <= 0 && dot(n, pts[s_j]) >= 0)
        swap(s_i, s_j);

    int long_edge_idx = 0;
    double max_len = sqr_length(pts[tri[1]] - pts[tri[0]]);
    for (int i = 1; i < 3; ++i)
    {
        double l = sqr_length(pts[tri[(i + 1) % 3]] - pts[tri[i]]);
        if (max_len < l)
        {
            max_len = l;
            long_edge_idx = i;
        }
    }

    int m = pts.size();
    pts.push_back(normalize(pts[tri[long_edge_idx]] + pts[tri[(long_edge_idx + 1) % 3]]));

    vector<Vec3i> triangles;
    for (int i = 0; i < 3; ++i)
    {

        if (long_edge_idx == i)
        {
            triangles.push_back(Vec3i(tri[i], m, s_i));
            triangles.push_back(Vec3i(m, tri[(i + 1) % 3], s_i));
            triangles.push_back(Vec3i(tri[(i + 1) % 3], m, s_j));
            triangles.push_back(Vec3i(m, tri[i], s_j));
        }
        else
        {
            triangles.push_back(Vec3i(tri[i], tri[(i + 1) % 3], s_i));
            triangles.push_back(Vec3i(tri[(i + 1) % 3], tri[i], s_j));
        }
    }
    return triangles;
}

void project_to_sphere(Manifold &m, const Vec3d &pn, double r)
{
    VertexAttributeVector<Vec3d> norms;
    for (auto v : m.vertices())
    {
        if (length(m.pos(v)) < 0.5)
            norms[v] = normal(m, v);
        else
            norms[v] = normalize(m.pos(v));
    }
    for (auto v : m.vertices())
    {
        m.pos(v) = pn + r * norms[v];
    }
};

void symmetrize_triangles(Manifold &m, VertexID v1, VertexID v2)
{
    for (auto h : m.halfedges())
    {
        auto w = m.walker(h);
        if (w.vertex() == v1 && w.opp().vertex() == v2)
        {
            m.split_edge(h);
            break;
        }
    }
    vector<FaceID> face_list;
    for (auto f : m.faces())
        face_list.push_back(f);
    for (FaceID f : face_list)
        m.split_face_by_vertex(f);
};

void symmetrize_tetrahedron(Manifold &m, VertexID v1, VertexID v2)
{
    for (auto h : m.halfedges())
    {
        auto w = m.walker(h);
        VertexID va = InvalidVertexID, vb = InvalidVertexID;
        if (w.vertex() == v1 && w.opp().vertex() == v2)
        {
            va = w.vertex();
            vb = w.opp().vertex();
        }
        else if (w.next().vertex() == v1 && w.opp().next().vertex() == v2)
        {
            va = w.next().vertex();
            vb = w.opp().next().vertex();
        }
        if (va == v1 && vb == v2)
        {
            VertexID v_new = m.split_edge(h);
            m.split_face_by_edge(w.face(), v_new, w.next().vertex());
            m.split_face_by_edge(w.opp().face(), v_new, w.opp().next().next().vertex());
        }
    }
};

void construct_bnps(HMesh::Manifold &m_out,
                    const Geometry::AMGraph3D &g,
                    Util::AttribVec<NodeID, FaceSet> &node2fs,
                    VertexAttributeVector<NodeID> &vertex2node,
                    const vector<double> &r_arr,
                    bool use_symmetry)
{

    map<int, pair<NodeID, NodeID>> spts2branch;
    map<int, VertexID> spts2vertexid;

    for (auto n : g.node_ids())
    {

        auto N = g.neighbors(n);
        if (N.size() > 2)
        {
            Manifold m;
            int node_vertex_count = 0;
            Vec3d pn = g.pos[n];

            vector<Vec3d> spts;
            int spts_vertex_count = 0;

            spts2branch.clear();
            spts2vertexid.clear();

            for (auto nn : N)
            {
                Vec3d pnn = g.pos[nn];
                spts.push_back(normalize(pnn - pn));

                auto spts_value = std::make_pair(n, nn);
                auto spts_key = spts_vertex_count;
                spts2branch.insert(std::make_pair(spts_key, spts_value));
                spts_vertex_count++;
            }

            // If we are supposed to symmetrize, we try to find symmetry pairs
            vector<pair<int, int>> npv;
            if (use_symmetry && N.size() < 6)
                npv = symmetry_pairs(g, n, 0.1);

            std::vector<CGLA::Vec3i> stris = SphereDelaunay(spts);

            if (N.size() == 5)
            {
                int n_ghosts = add_ghosts(stris, spts, 0.25);
                if (npv.size() > 0 && n_ghosts < 2)
                {
                    spts.resize(5);
                    stris = points_to_octahedron(spts, npv[0].first, npv[0].second);
                }
                else if (n_ghosts > 0)
                    stris = SphereDelaunay(spts);
            }
            else if (N.size() > 5 || npv.size() == 0)
            {
                int n_ghosts = add_ghosts(stris, spts, 0.75);
                if (n_ghosts > 0)
                    stris = SphereDelaunay(spts);
            }

            // Finally, we construct the BNP mesh from the triangle set.
            for (auto tri : stris)
            {
                vector<Vec3d> triangle_pts;
                for (int i = 0; i < 3; ++i)
                {
                    triangle_pts.push_back(spts[tri[i]]);
                    node_vertex_count++;
                }
                m.add_face(triangle_pts);
            }
            stitch_mesh(m, 1e-10);

            // Build mapping from spts to vertices in the BNP mesh
            for (auto v : m.vertices())
                for (int i = 0; i < spts.size(); i++)
                    if (sqr_length(m.pos(v) - spts[i]) < 0.0001)
                        spts2vertexid.insert(std::make_pair(i, v));

            // Refine the BNP mesh by splitting edges and faces if
            // symmetry pairs were found
            if (npv.size() > 0)
            {
                if (N.size() == 3)
                {
                    VertexID v1 = spts2vertexid[npv[0].first];
                    VertexID v2 = spts2vertexid[npv[0].second];
                    symmetrize_triangles(m, v1, v2);
                }
                else if (N.size() == 4)
                {
                    VertexID v1 = spts2vertexid[npv[0].first];
                    VertexID v2 = spts2vertexid[npv[0].second];
                    symmetrize_tetrahedron(m, v1, v2);
                }
            }

            // Project the BNP mesh to the sphere, make all vertices
            // valency 4, and do one step of catmull clark to make
            // the mesh a quadrilateral only mesh.
            project_to_sphere(m, pn, r_arr[n]);
            quad_valencify(m);
            id_preserving_cc(m);

            for (int i = 0; i < spts.size(); i++)
            {
                auto key = spts2branch.find(i)->second;
                auto value = m.pos(spts2vertexid.find(i)->second);
                branch2vert.insert(std::make_pair(key, value));
            }
            m.cleanup();

            size_t no_faces_before_merge = m_out.no_faces();
            size_t no_vertices_before_merge = m_out.allocated_vertices();

            m_out.merge(m);
            for (auto f : m_out.faces())
                if (f.index >= no_faces_before_merge)
                    node2fs[n].insert(f);
            for (auto v : m_out.vertices())
                if (v.index >= no_vertices_before_merge)
                    vertex2node[v] = n;
        }
    }
}

void merge_branch_faces(HMesh::Manifold &m, const Geometry::AMGraph3D &g,
                        const Util::AttribVec<NodeID, FaceSet> &node2fs)
{

    VertexID v;
    FaceID f, face_1, face_2;
    int branch_degree;
    HalfEdgeID boundary_edge_1;
    HalfEdgeID boundary_edge_2;

    for (auto n : g.node_ids())
    {
        auto N = g.neighbors(n);

        // for all branch nodes

        if (N.size() > 2)
        {

            // for each outgoing arc

            for (auto nn : N)
            {

                auto key = std::make_pair(n, nn);

                branch_degree = branchdeg.find(key)->second;

                f = branch_best_face.find(key)->second;

                v = branch_best_vertex.find(key)->second;

                if(valency(m,v) == branch_degree) {

                    HalfEdgeID ref_he;
                    VertexID ref_v;

                    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                        ref_he = w.halfedge();
                        if(m.walker(ref_he).vertex() == v)
                            ref_v = m.walker(ref_he).opp().vertex();
                        else if (m.walker(ref_he).opp().vertex() == v)
                            ref_v = m.walker(ref_he).vertex();
                    }

                    FaceID f = m.merge_one_ring(v);

                    if (m.in_use(ref_v))
                        one_ring_face_vertex[f] = ref_v;
                    else
                        one_ring_face_vertex[f] = InvalidVertexID;

                    branchface.insert(std::make_pair(key,f));
                    branch_best_vertex[key] = InvalidVertexID;
                    continue;

                }

                branchface.insert(std::make_pair(key,f));

                for(int i = 0; i < branch_degree - 1; i++) {

                    Walker w_f = m.walker(f);
                    Walker w_h = m.walker(w_f.halfedge());


                    HalfEdgeID w_start = w_f.halfedge();

                    do {

                        if(w_h.vertex() == v) {
                            boundary_edge_1 = w_h.halfedge();
                        }
                        if(w_h.opp().vertex() == v) {
                            boundary_edge_2 = w_h.halfedge();
                        }
                        w_h=w_h.next();

                    } while(w_h.halfedge() != w_start);

                    Walker b_e_1 = m.walker(boundary_edge_1);
                    Walker b_e_2 = m.walker(boundary_edge_2);
                    face_1 = b_e_1.opp().face();
                    face_2 = b_e_2.opp().face();


                    if( ! m.in_use(face_1)) {
                        m.merge_faces(f, boundary_edge_2);
                        continue;
                    }
                    if( ! m.in_use(face_2)) {
                        m.merge_faces(f, boundary_edge_1);
                        continue;
                    }


                    if(face_dist(m, g, n, nn, face_1) < face_dist(m, g, n, nn, face_2))
                        m.merge_faces(f, boundary_edge_1);
                    else
                        m.merge_faces(f, boundary_edge_2);

                }
            }
        }
    }

    return;
}

//Bridging Functions

vector<pair<VertexID, VertexID>> face_match_one_ring(const HMesh::Manifold& m, FaceID &f0, FaceID &f1,
                                                     const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {

    vector<pair<VertexID, VertexID> > connections;
    if(!m.in_use(f0) || !m.in_use(f1))
        return connections;

    VertexID face_vertex_0 = one_ring_face_vertex[f0];
    VertexID face_vertex_1 = one_ring_face_vertex[f1];

    bool fv_flag = false;

    if(face_vertex_0 != InvalidVertexID && face_vertex_1 != InvalidVertexID)
        fv_flag = true;

    int loop0_index = 0, loop1_index = 0;

    vector<VertexID> loop0;

    int count = 0;

    circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
        loop0.push_back(v);
        if(v == face_vertex_0)
            loop0_index = count;
        count++;
    }) );

    vector<VertexID> loop1;
    count = 0;

    circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
        loop1.push_back(v);
        if(v == face_vertex_1)
            loop1_index = count;
        count++;
    }) );

    size_t L0= loop0.size();
    size_t L1= loop1.size();

    if (L0 != L1)
        return connections;

    size_t L = L0;

    int j_off_min_len = -1;

    if(one_ring_face_vertex[f0] == InvalidVertexID || one_ring_face_vertex[f1] == InvalidVertexID) {

        float min_len = FLT_MAX;

        for(int j_off = 0; j_off < L; j_off = j_off + 1) {
            Vec3d bridge_edge_i, bridge_edge_j;
            float len = 0;

            for(int i=0;i<L;++i) {
                len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
            }
            if(len < min_len)   {
                j_off_min_len = j_off;
                min_len = len;
            }
        }

        int found_flag = 0;

        for (int i = 0; i < L; i++) {

            VertexID v0 = loop0[i];
            VertexID v1 = loop1[(L + j_off_min_len - i)%L];

            if(face_vertex_1 == v1) {
                Walker w = m.walker(f0);
                one_ring_face_vertex[f0] = v0;
                one_ring_face_vertex[m.walker(w.halfedge()).opp().face()] = v0;
                found_flag = 1;
            }

            else if (face_vertex_0 == v0) {
                Walker w = m.walker(f1);
                one_ring_face_vertex[f1] = v1;
                one_ring_face_vertex[m.walker(w.halfedge()).opp().face()] = v1;
                found_flag = 1;
            }
        }

        if(found_flag == 0 && (one_ring_face_vertex[f0] != InvalidVertexID || one_ring_face_vertex[f1] != InvalidVertexID)) {
            connections.clear();
            return connections;
        }
    }
    else {

        for(int j_off = 0; j_off < L; j_off = j_off + 1) {

            bool center_match = false;

            for(int i=0;i<L;++i) {
                if(loop0[i] == one_ring_face_vertex[f0] && loop1[(L + j_off - i)%L] == one_ring_face_vertex[f1])
                    center_match = true;

            }
            if(center_match) {
                j_off_min_len = j_off;
            }

        }
        float min_len = FLT_MAX;
        for(int j_off = j_off_min_len; j_off < 2*L; j_off = j_off + 2) {

            float len = 0;

            for(int i=0;i<L;++i)
                len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));


            if(len < min_len)   {
                j_off_min_len = j_off;
                min_len = len;
            }

        }
    }

    for(int i=0;i<L;++i)
        connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));


    return connections;
}

FaceID find_bridge_face(const HMesh::Manifold &m_out,
                        const Geometry::AMGraph3D& g, NodeID start_node, NodeID next_node, Util::AttribVec<NodeID, FaceSet>& node2fs) {

    FaceID f0 = InvalidFaceID;

    auto best_face = [&](NodeID n, NodeID nn) {
        Vec3d pn = g.pos[n];
        Vec3d pnn = g.pos[nn];
        Vec3d v_n_nn = pnn - pn;
        double d_max = -1000;
        FaceID f_max = InvalidFaceID;
        v_n_nn = normalize(v_n_nn);
        for(auto f: node2fs[n]) {
            double d = dot(v_n_nn, normal(m_out, f));
            if(d> d_max) {
                f_max = f;
                d_max = d;
            }

        }
        node2fs[n].erase(f_max);
        return f_max;
    };

    if(g.neighbors(start_node).size()>2)  {

        auto key = std::make_pair(start_node,next_node);
        f0 = branchface.find(key)->second;
        VertexID v0 = branch_best_vertex.find(key)->second;
    }
    else
        f0 = best_face(start_node,next_node);

    return f0;

}


void skeleton_aware_smoothing(const Geometry::AMGraph3D& g,
                              Manifold& m_out,
                              const VertexAttributeVector<NodeID>& vertex2node,
                              const vector<double>& node_radii) {
    const int N_iter = 50;
    for (int iter=0;iter<N_iter; ++iter) {


        Util::AttribVec<AMGraph::NodeID,Vec3d> barycenters(g.no_nodes(), Vec3d(0));
        Util::AttribVec<AMGraph::NodeID,int> cluster_cnt(g.no_nodes(), 0);
        for(auto v: m_out.vertices()) {
            NodeID n = vertex2node[v];
            barycenters[n] += m_out.pos(v);
            cluster_cnt[n] += 1;
        }

        for(auto n: g.node_ids()) 
            barycenters[n] /= cluster_cnt[n];

        auto new_pos = m_out.positions_attribute_vector();
        for (auto v: m_out.vertices()) {
            NodeID n = vertex2node[v];
            double w_sum = 0.0;
            Vec3d lap(0.0);
            Vec3d p0 = m_out.pos(v);
            for (auto vn: m_out.incident_vertices(v)) {
                double w = vertex2node[vn] == n ? 1 : 0.25;
                lap += w * (m_out.pos(vn)-p0);
                w_sum += w;
            }
            lap /= w_sum;
            Vec3d dir = cond_normalize(m_out.pos(v) + 0.5 * lap - barycenters[n]);
            Vec3d norm = normal(m_out, v);
            if (dot(dir, norm)<0.0)
                dir = norm;
            dir = cond_normalize(dir);
            double r = node_radii[n] * sqrt(g.valence(n)/2.0);
            new_pos[v] = 0.5 * (dir * r + g.pos[n] + m_out.pos(v));
        }
        m_out.positions_attribute_vector() = new_pos;
    }


}

void bridge_branch_node_meshes(Manifold& m_out, const AMGraph3D& g, Util::AttribVec<NodeID, FaceSet>& node2fs) {
    for(auto f_id: m_out.faces())
        if(one_ring_face_vertex.find(f_id) == one_ring_face_vertex.end())
            one_ring_face_vertex[f_id] = InvalidVertexID;

    for (auto n: g.node_ids()) {
        FaceID f0 = InvalidFaceID;
        VertexID v0, v1;

        auto N = g.neighbors(n);

        for(auto nn: N) {
            auto key = std::make_pair(n,nn);
            f0 = branchface.find(key)->second;

            NodeID start_node = n;
            NodeID next_node = nn;

            vector<NodeID> nbd_list = next_neighbours(g, start_node, next_node);

            do
            {
                FaceID f0 = find_bridge_face(m_out, g, start_node, next_node, node2fs);
                FaceID f1 = find_bridge_face(m_out, g, next_node, start_node, node2fs);

                nbd_list = next_neighbours(g, start_node, next_node);
                using VertexPair = pair<VertexID, VertexID>;
                vector<VertexPair> connections;

                if (not(f0 == InvalidFaceID || f1 == InvalidFaceID))
                    connections = face_match_one_ring(m_out, f0, f1, g, start_node, next_node);

                if (connections.size() != 0)
                    m_out.bridge_faces(f0, f1, connections);
                else
                    break;

                start_node = next_node;
                if (nbd_list.size() == 1)
                    next_node = nbd_list[0];

            } while (nbd_list.size() == 1);
        }
    }

}


//Main functions

HMesh::Manifold graph_to_FEQ(const Geometry::AMGraph3D& g, const vector<double>& _node_radii, bool use_symmetry) {

    Manifold m_out;
    Util::AttribVec<NodeID, FaceSet> node2fs;
    clear_global_arrays();
    double r = g.average_edge_length();
    vector<double> node_radii;
    node_radii.resize(g.no_nodes());
    for(auto n : g.node_ids()) {
        double l = r;
        for (auto m: g.neighbors(n))
            l = min(l, sqrt(g.sqr_dist(n, m)));
        node_radii[n] = 0.25*l;
    }

    VertexAttributeVector<NodeID> vertex2node(AMGraph::InvalidNodeID);
    construct_bnps(m_out, g, node2fs, vertex2node, node_radii, use_symmetry);
    init_branch_degree(m_out, g, node2fs);
    init_branch_face_pairs(m_out, g, node2fs);
    merge_branch_faces(m_out, g, node2fs);
    val2nodes_to_face_pairs(g, m_out, node2fs, vertex2node, node_radii);
    bridge_branch_node_meshes(m_out, g, node2fs);
    quad_mesh_leaves(m_out, vertex2node);
    skeleton_aware_smoothing(g, m_out, vertex2node, _node_radii);
    m_out.cleanup();
    return m_out;
}
