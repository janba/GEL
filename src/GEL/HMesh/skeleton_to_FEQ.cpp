#include <array>
#include <cmath>
#include <algorithm>
#include <unordered_map>

#include <GEL/CGLA/CGLA.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/SphereDelaunay.h>
#include <GEL/HMesh/quad_valencify.h>

using namespace Geometry;
using namespace CGLA;
using namespace HMesh;
using namespace std;
using NodeID = AMGraph::NodeID;


struct HashVec3d {
    size_t operator()(const CGLA::Vec3d& v) const noexcept {
        const size_t k = 83492791;
        size_t h = *reinterpret_cast<const size_t*>(&v[0]);
        for(unsigned int i=1; i<v.get_dim(); ++i) {
            h *= k;
            h += *reinterpret_cast<const size_t*>(&v[i]);
        }
        return h;
    }
};

// Initialize global arrays

struct BranchMeshInfo {
    HMesh::FaceID face;
    HMesh::VertexID vertex;
};

using BranchMeshMap = map<pair<NodeID,NodeID>,BranchMeshInfo>;
using Face2VertexMap = map<FaceID, VertexID>;

//Mesh util functions

void quad_mesh_leaves(HMesh::Manifold& m, VertexAttributeVector<NodeID>& vertex2node, Face2VertexMap& one_ring_face_vertex) {
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


//Functions for constructing / editing mesh elements from skeletal nodes

vector<FaceID> create_face_pair(Manifold& m, const Vec3d& pos, const Mat3x3d& _R, int axis,
                                Face2VertexMap& one_ring_face_vertex) {
    constexpr int num_sides = 8;
    Mat3x3d R = _R;
    double det = determinant(R);
    if(abs(det) > 1e-6)
        if(det<0) {
            Mat3x3d M = identity_Mat3x3d();
            M[2][2] = -1;
            R = R * M;
        }

    std::vector<FaceID> fvec;
    std::array<Vec3d, num_sides> pts;
    double angle = 0.0;
    for(int i = 0; i < num_sides; i++)
    {
        Vec3d _p = 0.5*Vec3d(0, cos(angle), sin(angle));
        angle += 2*M_PI/num_sides;
        Vec3d p(0);
        p[(0+axis)%3] += _p[0];
        p[(1+axis)%3] += _p[1];
        p[(2+axis)%3] += _p[2];
        pts[i] = R*p+pos;
    }
    fvec.push_back(m.add_face(pts));
    std::ranges::reverse(pts);
    fvec.push_back(m.add_face(pts));

    for (auto f: fvec) 
        one_ring_face_vertex[f] = m.walker(fvec[0]).vertex();

    return fvec;
}

void val2nodes_to_face_pairs(const Geometry::AMGraph3D &g, HMesh::Manifold &mani,
                             BranchMeshMap& branch2mesh_map,
                             VertexAttributeVector<NodeID> &vertex2node,
                             Face2VertexMap& one_ring_face_vertex,
                             const vector<double> &r)
{
    Util::AttribVec<NodeID, int> touched(g.no_nodes(), 0);
    Util::AttribVec<NodeID, Mat3x3d> warp_frame(g.no_nodes(), identity_Mat3x3d());
    Util::AttribVec<NodeID, Mat3x3d> mid_frame(g.no_nodes(), identity_Mat3x3d());
    Util::AttribVec<NodeID, int> max_idx(g.no_nodes(), 0);

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
                        size_t N_size = g.neighbors(m).size();
                        
                        Vec3d vect = normalize(g.pos[m] - g.pos[n]);
                        if (N_size == 2)
                        {
                            auto N_m = g.neighbors(m);
                            NodeID o = N_m[0]==n ? N_m[1] : N_m[0];
                            vect = normalize(vect + normalize(g.pos[o]-g.pos[m]));
                        }
                        Vec3d warp_v = warp_frame[n] * vect;

                        double max_sgn = sign(warp_v[0]);
                        double max_val = abs(warp_v[0]);
                        for (int i = 1; i < 3; ++i)
                        {
                            if (abs(warp_v[i]) > max_val)
                            {
                                max_sgn = sign(warp_v[i]);
                                max_val = abs(warp_v[i]);
                                max_idx[m] = i;
                            }
                        }
                        auto v_target = max_sgn * vect;
                        Quatd q;
                        q.make_rot((warp_frame[n])[max_idx[m]], v_target);
                        warp_frame[m] = transpose(q.get_Mat3x3d() * transpose(warp_frame[n]));

                    }
            }
        }
    
    for (auto m: g.node_ids())
    {
        size_t N_size = g.neighbors(m).size();
        if (N_size == 1 || N_size == 2)
        {
            Vec3d s(r[m]);
            Mat3x3d S = scaling_Mat3x3d(s);
            auto face_list = create_face_pair(mani, g.pos[m], transpose(warp_frame[m]) * S, max_idx[m], one_ring_face_vertex);
            stitch_mesh(mani, 1e-10);
            vector<Vec3d> norm;
            for (auto f : face_list)
            {
                for (auto v : mani.incident_vertices(f))
                    vertex2node[v] = m;
                norm.push_back(normal(mani, f));
            }
            int idx_taken=-1;
            for (auto mm: g.neighbors(m)) {
                Vec3d m_mm_v = cond_normalize(g.pos[mm] - g.pos[m]);
                int idx = (dot(norm[0], m_mm_v) > dot(norm[1], m_mm_v)) ? 0 : 1;
                if (idx==idx_taken)
                    idx = 1-idx;
                idx_taken = idx;
                branch2mesh_map[make_pair(m,mm)].face = face_list[idx];
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
        Vec3d n = cond_normalize(cross(v1, v2));
        ghost_pts.push_back(n);
    }

    /* Next, we cluster the ghost points. This is because in flatish
     configurations we could have several quite similar ghost points. */
    vector<int> cluster_id(ghost_pts.size(), -1);
    int max_id = 0;
    for (int i = 0; i < ghost_pts.size(); ++i)
    {
        if (cluster_id[i] == -1)
            cluster_id[i] = max_id++;
        for (int j = i + 1; j < ghost_pts.size(); ++j)
            if (cluster_id[j] == -1)
                if (dot(ghost_pts[i], ghost_pts[j]) > thresh)
                    cluster_id[j] = cluster_id[i];
    }

    vector<Vec3d> ghost_pts_new(max_id, Vec3d(0));
    for (int i = 0; i < ghost_pts.size(); ++i)
        ghost_pts_new[cluster_id[i]] += ghost_pts[i];

    /* Finally, we cull ghost points too close to an existing non-ghost point. */
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


vector<Vec3i> five_points_to_octahedron(vector<Vec3d> &pts, int s_i, int s_j)
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

tuple<Manifold, BranchMeshMap, VertexAttributeVector<NodeID>>
construct_bnps(const Geometry::AMGraph3D &g,
               const vector<double> &r_arr,
               bool use_symmetry)
{
    Manifold m_out;
    BranchMeshMap branch2mesh_map;
    VertexAttributeVector<NodeID> vertex2node(AMGraph::InvalidNodeID);

    unordered_map<Vec3d, pair<NodeID,NodeID>,HashVec3d> pos_to_branch;
    for (auto n : g.node_ids())
    {
        auto N = g.neighbors(n);
        if (N.size() > 2)
        {
            Manifold m;
            Vec3d pn = g.pos[n];

            vector<Vec3d> spts;
            vector<pair<NodeID, NodeID>> spts2branch;
            
            for (auto nn : N) {
                Vec3d pnn = g.pos[nn];
                spts.push_back(normalize(pnn - pn));
                spts2branch.push_back(std::make_pair(n, nn));
            }

            // If we are supposed to symmetrize, we try to find symmetry pairs
            vector<pair<int, int>> npv;
            if (N.size() < 6)
                npv = symmetry_pairs(g, n, 0.1, !use_symmetry);


            std::vector<CGLA::Vec3i> stris = SphereDelaunay(spts);

            if (N.size() < 5) {
                if (npv.size() == 0) {
                    int n_ghosts = add_ghosts(stris, spts, 0.8);
                    if (n_ghosts > 0)
                        stris = SphereDelaunay(spts);
                    else cout << "NPV size " << npv.size() << " " << N.size() <<  " " << n_ghosts << endl;
                }
            }
            else if (N.size() == 5) {
                int n_ghosts = add_ghosts(stris, spts, 0.25);
                if (npv.size() > 0 && n_ghosts < 2)
                {
                    spts.resize(5);
                    stris = five_points_to_octahedron(spts, npv[0].first, npv[0].second);
                }
                else if (n_ghosts > 0)
                    stris = SphereDelaunay(spts);
            }
            else if (N.size() > 5) {
                int n_ghosts = add_ghosts(stris, spts, 0.8);
                if (n_ghosts > 0)
                    stris = SphereDelaunay(spts);
            }

            // Finally, we construct the BNP mesh from the triangle set.
            VertexAttributeVector<int> vertexid2spts;
            for (auto tri : stris)
            {
                std::array<Vec3d, 3> triangle_pts;
                for (int i = 0; i < triangle_pts.size(); ++i)
                {
                    triangle_pts[i] = spts[tri[i]];
                }
                FaceID f = m.add_face(triangle_pts);
                int i=0;
                for (auto v: m.incident_vertices(f))
                    vertexid2spts[v] = tri[i++];
            }

            stitch_mesh(m, 1e-10);

            vector<VertexID> spts2vertexid(spts.size());
            // Build mapping from spts to vertices in the BNP mesh
            for (auto v : m.vertices()) {
                int i = vertexid2spts[v];
                spts2vertexid[i] = v;
            }

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
            
            for(auto v: m.vertices()) {
                if (valency(m, v) != 4) {
                    cout << "bad bad face: " << endl;
                    cout << "stats: " << m.no_vertices() << ", " << m.no_faces() << endl;
                    cout << "incident edges: " << N.size() << endl;
                    cout << "symmetries: " << npv.size() << endl;
                    cout << "valencies: ";
                    for (auto vv: m.vertices())
                        cout << valency(m, vv) <<  " ";
                    cout << endl;
                    break;
                }
            }
            
            cc_split(m);

            for (int i = 0; i < spts.size(); i++)
                pos_to_branch[m.pos(spts2vertexid[i])] = spts2branch[i];

            m.cleanup();

            size_t no_faces_before_merge = m_out.no_faces();
            size_t no_vertices_before_merge = m_out.allocated_vertices();

            m_out.merge(m);

            for (auto v : m_out.vertices())
                if (v.index >= no_vertices_before_merge)
                    vertex2node[v] = n;
        }
    }
    for(auto v: m_out.vertices())
        branch2mesh_map[pos_to_branch[m_out.pos(v)]].vertex = v;

    return { m_out, branch2mesh_map, vertex2node };
}

Face2VertexMap merge_branch_faces(const Geometry::AMGraph3D &g,
                        HMesh::Manifold &m,
                        BranchMeshMap& branch2mesh_map)
{
    Face2VertexMap one_ring_face_vertex;
    for (auto n : g.node_ids())
    {
        auto N = g.neighbors(n);

        // for all branch nodes

        if (N.size() > 2)
        {
            // for each outgoing arc
            for (auto nn : N)
            {
                auto& b2mm = branch2mesh_map[make_pair(n,nn)];
                VertexID v = b2mm.vertex;

                if(valency(m,v) == 4) {

                    HalfEdgeID ref_he;
                    VertexID ref_v;

                    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                        ref_he = w.halfedge();
                        if(m.walker(ref_he).vertex() == v)
                            ref_v = m.walker(ref_he).opp().vertex();
                        else if (m.walker(ref_he).opp().vertex() == v)
                            ref_v = m.walker(ref_he).vertex();
                    }

                    FaceID f_to_merge = m.merge_one_ring(v);

                    if (m.in_use(ref_v))
                        one_ring_face_vertex[f_to_merge] = ref_v;
                    else
                        one_ring_face_vertex[f_to_merge] = InvalidVertexID;

                    b2mm.face = f_to_merge;
                }
            }
        }
    }

    return one_ring_face_vertex;
}

//Bridging Functions

vector<pair<VertexID, VertexID>> face_match_one_ring(const HMesh::Manifold& m, FaceID &f0, FaceID &f1, Face2VertexMap& one_ring_face_vertex) {
    
    vector<pair<VertexID, VertexID> > connections;
    if(!m.in_use(f0) || !m.in_use(f1)) {
        cout << "one face unused" << endl;
        return connections;
    }
    
    VertexID face_vertex_0 = one_ring_face_vertex[f0];
    VertexID face_vertex_1 = one_ring_face_vertex[f1];
    
    int loop0_index = 0, loop1_index = 0;
    
    vector<VertexID> loop0;
    
    int count = 0;
    
    circulate_face_ccw(m, f0, [&](VertexID v) {
        loop0.push_back(v);
        if (v == face_vertex_0) loop0_index = count;
        count++;
    });
    
    vector<VertexID> loop1;
    count = 0;
    
    circulate_face_ccw(m, f1, [&](VertexID v) {
        loop1.push_back(v);
        if (v == face_vertex_1) loop1_index = count;
        count++;
    });
    
    size_t L0= loop0.size();
    size_t L1= loop1.size();
    
    if (L0 != L1) {
        cout << "loop sizes " << L0 << " " << L1 << endl;   
        return connections;
    }
    
    size_t L = L0;
    
    int j_off_min_len = -1;
    
    for(int j_off = 0; j_off < L; j_off = j_off + 1) {
        bool center_match = false;
        
        for(int i=0;i<L;++i) {
            if(loop0[i] == one_ring_face_vertex[f0] && loop1[(L + j_off - i)%L] == one_ring_face_vertex[f1])
                center_match = true;
        }
        if(center_match)
            j_off_min_len = j_off;
    }
    double min_len = FLT_MAX;
    for(int j_off = j_off_min_len; j_off < 2*L; j_off = j_off + 2) {
        double len = 0;
        for(int i=0;i<L;++i)
            len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
        
        if(len < min_len)   {
            j_off_min_len = j_off;
            min_len = len;
        }
        
    }
    
    for(int i=0;i<L;++i)
        connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
    
    return connections;
}



void bridge_branch_node_meshes(const AMGraph3D& g,
                               Manifold& m_out,
                               BranchMeshMap& branch2mesh_map,
                               Face2VertexMap& one_ring_face_vertex) {
    for(auto f_id: m_out.faces())
        if(one_ring_face_vertex.find(f_id) == one_ring_face_vertex.end())
            one_ring_face_vertex[f_id] = InvalidVertexID;

    Util::AttribVec<NodeID,int> visited(0);
    for (auto n: g.node_ids()) {
        visited[n] = 1;

        for(auto nn:  g.neighbors(n)) 
            if (!visited[nn]) {
                FaceID f0 = branch2mesh_map[make_pair(n,nn)].face;
                FaceID f1 = branch2mesh_map[make_pair(nn,n)].face;

                using VertexPair = pair<VertexID, VertexID>;
                vector<VertexPair> connections;

                if (not(f0 == InvalidFaceID || f1 == InvalidFaceID))
                    connections = face_match_one_ring(m_out, f0, f1, one_ring_face_vertex);

                if (connections.size() != 0)
                    m_out.bridge_faces(f0, f1, connections);
                else
                    cout << "FAILED TO CONNECT " << f0.get_index() << " " << f1.get_index() << endl;
        }
    }

}

void skeleton_aware_smoothing(const Geometry::AMGraph3D& g,
                              Manifold& m_out,
                              const VertexAttributeVector<NodeID>& vertex2node,
                              const vector<double>& node_radii) {
    const int N_dir_idx = 10;
    for (int dir_idx=0;dir_idx<N_dir_idx; ++dir_idx) {
        auto new_pos = m_out.positions_attribute_vector();
        for (auto v: m_out.vertices()) {
            NodeID n = vertex2node[v];
            double w_sum = 0.0;
            Vec3d lap(0.0);
            Vec3d p0 = m_out.pos(v);
            for (auto vn: m_out.incident_vertices(v)) {
                lap += m_out.pos(vn)-p0;
                w_sum += 1.0;
            }
            lap /= w_sum;
            new_pos[v] = m_out.pos(v) + lap;
        }
        m_out.positions_attribute_vector() = new_pos;

        Util::AttribVec<AMGraph::NodeID,Vec3d> barycenters(g.no_nodes(), Vec3d(0));
        Util::AttribVec<AMGraph::NodeID,int> cluster_cnt(g.no_nodes(), 0);
        for(auto v: m_out.vertices()) {
            NodeID n = vertex2node[v];
            barycenters[n] += m_out.pos(v);
            cluster_cnt[n] += 1;
        }
        for(auto n: g.node_ids())
            barycenters[n] /= cluster_cnt[n];

        for (auto v: m_out.vertices()) {
            NodeID n = vertex2node[v];
            Vec3d norm = normal(m_out, v);
            Vec3d dir = cond_normalize(m_out.pos(v) - barycenters[n]);
            double r = node_radii[n] * sqrt(g.valence(n)/2.0);
            // new_pos[v] = dir * r + g.pos[n];
            new_pos[v] = 0.7 * ((0.7*dir+0.3*norm) * r + g.pos[n]) + 0.3 * m_out.pos(v);
        }
        m_out.positions_attribute_vector() = new_pos;
    }


}



//Main functions



HMesh::Manifold graph_to_FEQ(const Geometry::AMGraph3D& g, const vector<double>& _node_radii, bool use_symmetry) {
    double r = g.average_edge_length();
    vector<double> node_radii;
    node_radii.resize(g.no_nodes());
    for(auto n : g.node_ids()) {
        double l = r;
        for (auto m: g.neighbors(n))
            l = min(l, sqrt(g.sqr_dist(n, m)));
        node_radii[n] = 0.25*l;
    }
    auto [m_out, branch2mesh_map, vertex2node] = construct_bnps(g, node_radii, use_symmetry);
    auto face_vertex = merge_branch_faces(g, m_out, branch2mesh_map);
    val2nodes_to_face_pairs(g, m_out, branch2mesh_map, vertex2node, face_vertex, node_radii);
    bridge_branch_node_meshes(g, m_out, branch2mesh_map, face_vertex);
    quad_mesh_leaves(m_out, vertex2node, face_vertex);
    skeleton_aware_smoothing(g, m_out, vertex2node, _node_radii);
    m_out.cleanup();
    return m_out;
}

std::vector<Vec3d> generateVectors(int N) {
    std::vector<Vec3d> vectors;
    auto rand_num = []() { return static_cast<double>(rand())/static_cast<double>(RAND_MAX); };

    for (int i = 0; i < N; ++i) {
        double inclination = acos(1 - 2 * rand_num());
        double azimuth = 2 * M_PI * rand_num();

        double x = sin(inclination) * cos(azimuth);
        double y = sin(inclination) * sin(azimuth);
        double z = cos(inclination);

        vectors.push_back(Vec3d(x, y, z));
    }

    return vectors;
}




void non_rigid_registration(HMesh::Manifold& m, const HMesh::Manifold& m_ref) {
    srand(0);

    auto rand_num = []() { return static_cast<double>(rand())/static_cast<double>(RAND_MAX); };
    const int N_dir = 47; // Number of directions
    // Generate random vectors distributed over the sphere.
    vector<Vec3d> dir = generateVectors(N_dir);
    static int call_no = 0;

    cout << __FILE__ << __LINE__ << " : " << call_no << endl;
    const int N_pts = m_ref.no_vertices(); // Number of points

    vector<Vec3d> m_pts;      // Points on manifold we deform
    vector<Vec3d> m_ref_pts;  // Target points on reference manifold
    FaceAttributeVector<vector<int>> m_pts_idx; // Point indices stored per face
    cout << __FILE__ << __LINE__ << " : " << call_no << endl;

    // Generate points on the reference manifolds. These are simply
    // its vertices.
    for(auto p: m_ref.vertices())
        m_ref_pts.push_back(m_ref.pos(p));
    cout << __FILE__ << __LINE__ << " : " << call_no << endl;

    // Now compute the total area of the deformable manifold, m, and 
    // also stored the areas of its faces in an attribute vector.
    double total_area = 0;
    FaceAttributeVector<double> face_area;
    for (const FaceID f: m.faces()) {
        double a = area(m, f);
        if (std::isnan(a))
            a = 0;
        face_area[f] = a;
        total_area += a;
    }
    // cout << __FILE__ << __LINE__ << " : " << call_no << endl;

    // Generate N_pts points randomly distributed over m.
    // The face areas are used to ensure an even distribution.
    for (int i=0;i<N_pts;++i) {
        FaceID f = InvalidFaceID;
        // cout << __FILE__ << __LINE__ << " : " << call_no << " " << f.get_index() << " "  << i << endl;

        while (f == InvalidFaceID) {
            double r = rand_num() * total_area;
            // cout << __FILE__ << __LINE__ << " : " << call_no << " " << f.get_index() << " "  << r << " " << total_area << endl;
            double a = 0;
            for (const FaceID ff: m.faces()) {
                a += face_area[ff];
                if (a+1e-6>r) {
                    f = ff;
                    break;
                }
            }
        }
        double w_sum = 0;
        Vec3d p(0);
        // The point is just a random affine combination of the polygon's
        // vertices.
        for (VertexID v: m.incident_vertices(f)) {
            double w = rand_num();
            p += w * m.pos(v);
            w_sum += w;
        }
        m_pts_idx[f].push_back(m_pts.size());
        m_pts.push_back(p/w_sum);
    }
    cout << __FILE__ << __LINE__ << " : " << call_no << endl;

    VertexAttributeVector<Vec3d> new_pos(Vec3d(0));
    for (int iter = 0; iter < N_dir; ++iter)
    {
        Vec3d V = dir[iter];

        vector<pair<double, int>> pts_m_1d;
        vector<pair<double, int>> pts_ref_1d;
        for (int j = 0; j < N_pts; ++j)
        {
            double d = dot(V, m_pts[j]);
            pts_m_1d.push_back(pair<double, int>(d, j));
            d = dot(V, m_ref_pts[j]);
            pts_ref_1d.push_back(pair<double, int>(d, j));
        }

        sort(pts_m_1d.begin(), pts_m_1d.end());
        sort(pts_ref_1d.begin(), pts_ref_1d.end());

        vector<Vec3d> displace(N_pts, Vec3d(0));
        for (int j = 0; j < N_pts; ++j)
        {
            int i_m = pts_m_1d[j].second;
            int i_ref = pts_ref_1d[j].second;
            displace[i_m] = m_ref_pts[i_ref]-m_pts[i_m];
//           displace[i_m] = V * dot(V, m_ref_pts[i_ref]-m_pts[i_m]);
        }

        for (VertexID v : m.vertices())
        {
            Vec3d p(0);
            int n = 0;
            for (FaceID f : m.incident_faces(v))
            {
                for (int i : m_pts_idx[f])
                {
                    p += displace[i];
                    ++n;
                }
            }
            new_pos[v] += p / n + m.pos(v);
        }
    }
    cout << __FILE__ << __LINE__ << " : " << call_no << endl;

    for (auto v: m.vertices())
        m.pos(v) = new_pos[v]/N_dir;
    cout << __FILE__ << __LINE__ << " : " << call_no << endl;
    call_no++;

}
