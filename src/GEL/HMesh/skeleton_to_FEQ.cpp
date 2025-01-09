#include <array>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#include <GEL/CGLA/CGLA.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/SphereDelaunay.h>
#include <GEL/HMesh/quad_valencify.h>
#include <GEL/HMesh/face_loop.h>

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
//    HMesh::VertexID vertex;
};

using BranchMeshMap = map<pair<NodeID,NodeID>,BranchMeshInfo>;
using Face2VertexMap = map<FaceID, VertexID>;
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
    vector<HalfEdgeID> base_edges;
    vector<HalfEdgeID> new_edges;

    for(auto f: m_in.faces())
        base_faces.push_back(f);

    for(auto h: m_in.halfedges())
        if (h < m_in.walker(h).opp().halfedge())
            base_edges.push_back(h);

    for(auto f: base_faces) {
        VertexID center_v = m_in.split_face_by_vertex(f);
        for (auto h: m_in.incident_halfedges(center_v))
            new_edges.push_back(h);
    }

    for (auto h: base_edges) {
        Walker w = m_in.walker(h);
        VertexID v1 = w.next().vertex();
        FaceID f1 = w.face();
        VertexID v2 = w.opp().next().vertex();
        FaceID f2 = w.opp().face();
        VertexID v = m_in.split_edge(h);
        m_in.split_face_by_edge(f1, v1, v);
        m_in.split_face_by_edge(f2, v2, v);
    }

    for (auto h_dissolve: new_edges)
        if(m_in.in_use(h_dissolve))
            m_in.merge_faces(m_in.walker(h_dissolve).face(), h_dissolve);
    return;

}

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

vector<FaceID> create_face_pair(Manifold& m, const Vec3d& pos, const Mat3x3d& R) {
    int num_sides = 4;

    vector<FaceID> fvec;
    vector<Vec3d> pts;
    for(int i = 0; i < num_sides; i++)
    {
        double angle = 2.0*i*M_PI/num_sides;
        Vec3d p = 0.5*Vec3d(cos(angle), sin(angle), 0);
        pts.push_back(R*p+pos);
    }
    fvec.push_back(m.add_face(pts));
    reverse(begin(pts), end(pts));
    fvec.push_back(m.add_face(pts));
    return fvec;
}

void val2nodes_to_face_pairs(const Geometry::AMGraph3D &g, HMesh::Manifold &mani,
                             BranchMeshMap& branch2mesh_map,
                             VertexAttributeVector<NodeID> &vertex2node,
                             const vector<double> &r)
{

    for (auto m: g.node_ids())
    {
        auto N = g.neighbors(m);
        size_t N_size = g.neighbors(m).size();
        if (N_size == 1 || N_size == 2)
        {
            double sgn = 1;
            Vec3d Z(0);
            for (auto n: N) {
                Z += sgn*(g.pos[n]-g.pos[m]);
                sgn = -sgn;
            }
            Z = normalize(Z);
            Vec3d X, Y;
            orthogonal(Z, X, Y);
            
            Mat3x3d M = r[m]*transpose(Mat3x3d(X, Y, Z));

            auto face_list = create_face_pair(mani, g.pos[m], M);
            stitch_mesh(mani, 1e-10);
            for (auto f : face_list)
            {
                for (auto v : mani.incident_vertices(f))
                    vertex2node[v] = m;
            }
            int idx_sum=0;
            for (auto mm: g.neighbors(m)) {
                Vec3d m_mm_v = g.pos[mm] - g.pos[m];
                Vec3d n0 = normal(mani, face_list[0]);
                Vec3d n1 = normal(mani, face_list[1]);
                
                int idx = (dot(n0, m_mm_v) > dot(n1, m_mm_v)) ? 0 : 1;
                branch2mesh_map[make_pair(m,mm)].face = face_list[idx];
            }
            if(idx_sum == 1)
                cout << "error: doubly assigned face in bridge node " << endl;
        }
    }

}
/**
 void val2nodes_to_face_pairs(const Geometry::AMGraph3D &g, HMesh::Manifold &mani,
                              BranchMeshMap& branch2mesh_map,
                              VertexAttributeVector<NodeID> &vertex2node,
                              const vector<double> &r)
 {
     Util::AttribVec<NodeID, int> touched(g.no_nodes(), 0);
     Util::AttribVec<NodeID, Mat3x3d> warp_frame(g.no_nodes(), identity_Mat3x3d());
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
                             NodeID o = next_neighbours(g, n, m)[0];
                             vect = normalize((g.pos[o]-g.pos[n]));
                         }
                         vect = normalize(vect);
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
             auto face_list = create_face_pair(mani, g.pos[m], transpose(warp_frame[m]) * S, max_idx[m]);
             stitch_mesh(mani, 1e-10);
             for (auto f : face_list)
             {
                 for (auto v : mani.incident_vertices(f))
                     vertex2node[v] = m;
             }
             int idx_sum=0;
             for (auto mm: g.neighbors(m)) {
                 Vec3d m_mm_v = g.pos[mm] - g.pos[m];
                 Vec3d n0 = normal(mani, face_list[0]);
                 Vec3d n1 = normal(mani, face_list[1]);
                 
                 int idx = (dot(n0, m_mm_v) > dot(n1, m_mm_v)) ? 0 : 1;
                 branch2mesh_map[make_pair(m,mm)].face = face_list[idx];
             }
             if(idx_sum == 1)
                 cout << "error: doubly assigned face in bridge node " << endl;
         }
     }

 }
**/

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
            
            for (auto nn : N) {
                Vec3d pnn = g.pos[nn];
                spts.push_back(normalize(pnn - pn));
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
                vector<Vec3d> triangle_pts;
                for (int i = 0; i < 3; ++i)
                {
                    triangle_pts.push_back(spts[tri[i]]);
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
            
            // Project the BNP mesh to the sphere and make all vertices
            // valency 4
            VertexAttributeVector<Vec3d> norms;
            for (auto v : m.vertices())
            {
                if (length(m.pos(v)) < 0.5)
                    norms[v] = normal(m, v);
                else
                    norms[v] = normalize(m.pos(v));
            }
            for (auto v : m.vertices())
                m.pos(v) = norms[v];
            
//            string path = "/Users/janba/GEL/src/demo/FEQ-Remeshing/";
//            time_t x = time(0) + random();
//            string file_name = path + "BNP_" + to_string(x) + ".obj";
//            obj_save(file_name, m);
//            HMesh::VertexAttributeVector<CGLA::Vec3f> vcol;
//            HMesh::HalfEdgeAttributeVector<CGLA::Vec3f> hcol;
//            HMesh::FaceAttributeVector<CGLA::Vec3f> fcol;
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
            
            Manifold m_dual;
            FaceAttributeVector<Vec3d> dual_verts;
            for (auto f: m.faces())
                dual_verts[f] = normalize(centre(m,f))*r_arr[n]+pn;
            for (auto v: m.vertices()) {
                vector<Vec3d> barycenters;
                for (auto f: m.incident_faces(v))
                    barycenters.push_back(dual_verts[f]);
                m_dual.add_face(barycenters);
            }
            stitch_mesh(m_dual, 1e-10);
//            while(collapse_double_crossed_quads(m_dual));
//            int cnt = split_double_crossed_quads(m_dual) + split_double_crossed_quads(m_dual);
//            cout << "Split " << cnt << " double crossed quad" << endl;
            m_dual.cleanup();
            
            size_t no_vertices_before_merge = m_out.allocated_vertices();
            size_t no_faces_before_merge = m_out.allocated_faces();
            m_out.merge(m_dual);

            for (auto v : m_out.vertices())
                if (v.index >= no_vertices_before_merge)
                    vertex2node[v] = n;
            
            FaceSet fset;
            for (auto f: m_out.faces())
                if (f.index >= no_faces_before_merge)
                    fset.insert(f);

            for (auto nn : N) {
                double max_dot=-1;
                FaceID max_dot_f = InvalidFaceID;
                Vec3d branch_dir = normalize(g.pos[nn] - pn);
                for (auto f: fset) {
                    Vec3d norm = normal(m_out, f);
                    double d = dot(branch_dir, norm);
                    if (d > max_dot) {
                        max_dot = d;
                        max_dot_f = f;
                    }
                }
                branch2mesh_map[std::make_pair(n, nn)].face = max_dot_f;
                fset.erase(max_dot_f);
            }
                    
                
        }
    }
    return { m_out, branch2mesh_map, vertex2node };
}

//Bridging Functions

vector<pair<VertexID, VertexID>> face_match_one_ring(const HMesh::Manifold& m, FaceID &f0, FaceID &f1) {
    
    vector<pair<VertexID, VertexID> > connections;
    if(!m.in_use(f0) || !m.in_use(f1))
        return connections;
    
    vector<VertexID> loop0;
    Vec3d n0 = normal(m, f0);
    circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
        loop0.push_back(v);
    }) );
    
    vector<VertexID> loop1;
    Vec3d n1 = normal(m, f1);
    circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
        loop1.push_back(v);
    }) );
    
    size_t L0= loop0.size();
    size_t L1= loop1.size();
    
    if (L0 != L1)
        return connections;
    
    size_t L = L0;
    
    int j_off_min_len = 0;
    double min_len = FLT_MAX;
    for(int j_off = 0; j_off < L; j_off += 1) {
        double len = 0;
        for(int i=0;i<L;++i) {
            Vec3d v = m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]);
            len += sqr_length(v-n0*dot(v,n0)) + sqr_length(v-n1*dot(v,n1));
        }
        if(len < min_len)   {
            j_off_min_len = j_off;
            min_len = len;
        }
        
    }
    
    for(int i=0;i<L;++i)
        connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
    
    return connections;
}

void align_branch_node_meshes(const AMGraph3D& g,
                               Manifold& m_out,
                               BranchMeshMap& branch2mesh_map) {
    BreadthFirstSearch bfs(g);
    for (auto n: g.node_ids()) {
        auto N = g.neighbors(n);
        if (N.size() >= 3)
            bfs.add_init_node(n);
    }
    while(bfs.Dijkstra_step());
    
    vector<pair<double, NodeID>> pri_node_vec;
    double d0 = 1.0;
    for (auto n: g.node_ids()) {
        if (g.neighbors(n).size()<3)
        {
            if(bfs.dist[n]>1e100) {
                bfs.dist[n] = d0;
                d0 = d0 + 1;
            }
            pri_node_vec.push_back(make_pair(bfs.dist[n], n));
        }
    }
    
    sort(begin(pri_node_vec), end(pri_node_vec));
    
    
    for (int iter=0;iter<5;++iter) {
        double w = 1.0;//0.2 + 0.5 * (iter/50.0);
        for (auto [priority, n]: pri_node_vec) {
            auto N = g.neighbors(n);
            for(auto nn: N)
                if (bfs.dist[nn] <= priority) {
//                    cout << "aligning " << n << " ( " << bfs.dist[n] << " ) and " << nn << " ( " << bfs.dist[nn] << " ) " << endl;
                    FaceID f0 = branch2mesh_map[make_pair(n,nn)].face;
                    FaceID f1 = branch2mesh_map[make_pair(nn,n)].face;
                    if (not(f0 == InvalidFaceID || f1 == InvalidFaceID)) {
                        auto connections = face_match_one_ring(m_out, f0, f1);
                        if (connections.size() != 0)
                        {
                            Vec3d c0 = barycenter(m_out, f0);
                            Vec3d Z = normalize(g.pos[nn]-g.pos[n]);
                            double alpha_sum = 0;
                            for (auto [v0, v1]: connections) {
                                Vec3d X = m_out.pos(v0) - c0;
                                Vec3d Y = cross(Z,X);
                                Vec3d p1 = m_out.pos(v1) - c0;
                                alpha_sum += atan2(dot(p1,Y), dot(p1,X));
                            }
                            double alpha = w * alpha_sum / 4;
                            for (auto v: m_out.incident_vertices(f0)) {
                                Vec3d X = m_out.pos(v) - c0;
                                Vec3d Y = cross(Z,X);
                                m_out.pos(v) = X * cos(alpha) + Y * sin(alpha) + c0;
                            }
                        }
                    }
                }
//                else {
//                    cout << "not aligning " << n << " ( " << bfs.dist[n] << " ) and " << nn << " ( " << bfs.dist[nn] << " ) " << endl;
//                }
        }
    }
}


void bridge_branch_node_meshes(const AMGraph3D& g,
                               Manifold& m_out,
                               BranchMeshMap& branch2mesh_map) {
    for (auto n: g.node_ids()) {
        auto N = g.neighbors(n);

        for(auto nn: N) {
            auto key = std::make_pair(n,nn);
            
            FaceID f0 = branch2mesh_map[make_pair(n,nn)].face;
            FaceID f1 = branch2mesh_map[make_pair(nn,n)].face;
            if (not(f0 == InvalidFaceID || f1 == InvalidFaceID)) {
                auto connections = face_match_one_ring(m_out, f0, f1);
                if (connections.size() != 0)
                    m_out.bridge_faces(f0, f1, connections);
            
            }
        }
    }
}

void skeleton_aware_smoothing(const Geometry::AMGraph3D& g,
                              Manifold& m_out,
                              const VertexAttributeVector<NodeID>& vertex2node,
                              const vector<double>& node_radii) {
    const int N_dir_idx = 5;
    for (int dir_idx=0; dir_idx<N_dir_idx; ++dir_idx) {
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
        m_out.positions = new_pos;
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
    val2nodes_to_face_pairs(g, m_out, branch2mesh_map, vertex2node, node_radii);
    align_branch_node_meshes(g, m_out, branch2mesh_map);
    bridge_branch_node_meshes(g, m_out, branch2mesh_map);
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
