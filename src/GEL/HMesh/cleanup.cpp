/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <map>
#include "cleanup.h"

#include "../CGLA/Vec3f.h"
#include "../CGLA/Vec3d.h"
#include "../Geometry/QEM.h"
#include "../Geometry/KDTree.h"

#include "refine_edges.h"
#include "Manifold.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    using namespace Geometry;
    
    namespace
    {
        /// get the minimum angle between 3 vertices
        float min_angle(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2)
        {
            Vec3d a = normalize(v1-v0);
            Vec3d b = normalize(v2-v1);
            Vec3d c = normalize(v0-v2);
            
            return min(dot(a, -c), min(dot(b, -a), dot(c, -b)));
        }
        
        /// need description
        float edge_error(const Manifold& m, HalfEdgeID h, const Vec3d& pa, const Vec3d& pb)
        {
            QEM q;
            Walker j = m.walker(h);
            
            FaceID f = j.face();
            if(f != InvalidFaceID)
                q += QEM(Vec3d(0), Vec3d(normal(m, f)));
            
            f = j.opp().face();
            if(f != InvalidFaceID)
                q += QEM(Vec3d(0), Vec3d(normal(m, f)));
            
            return q.error(pb - pa);
        }
        
        /// need description
        float vertex_error(const Manifold& m, VertexID v, const Vec3d& pb)
        {
            QEM q;
            Vec3d pa(m.pos(v));
            
            for(Walker vj = m.walker(v); !vj.full_circle(); vj = vj.circulate_vertex_cw()){
                FaceID f = vj.face();
                if(f != InvalidFaceID)
                    q += QEM(Vec3d(0), Vec3d(normal(m, f)));
            }
            return q.error(pb - pa);
        }
    }
    
    void remove_caps(Manifold& m, float thresh)
    {
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            Vec3d p[3];
            HalfEdgeID he[3];
            VertexID vh[3];
            
            // store ids of vertices and halfedges and vertex positions of face
            size_t n = 0;
            for(Walker fj = m.walker(*f); !fj.full_circle(); fj = fj.circulate_face_cw(), ++n){
                vh[n] = fj.vertex();
                p[n] = Vec3d(m.pos(vh[n]));
                
                // original face circulator implementation returned next halfedge, jumper doesn't. Can this be optimized?
                he[n] = fj.halfedge();
            }
            assert(n == 3);
            
            // calculate the edge lengths of face
            bool is_collapsed = false;
            Vec3d edges[3];
            for(size_t i = 0; i < 3; ++i){
                edges[i] = p[(i+1)%3] - p[i];
                float l = length(edges[i]);
                if(l < 1e-20)
                    is_collapsed = true;
                else
                    edges[i] /= l;
                
            }
            // an edge length was close to 1e-20, thus collapsed
            if(is_collapsed)
                continue;
            
            for(size_t i = 0; i < 3; ++i){
                float ang = acos(max(-1.0, min(1.0, dot(-edges[(i+2)%3], edges[i]))));
                
                // flip long edge of current face if angle exceeds cap threshold and result is better than current cap
                if(ang > thresh){
                    size_t iplus1 = (i+1)%3;
                    Vec3d edge_dir = edges[iplus1];
                    
                    Walker j = m.walker(he[iplus1]);
                    Vec3d v0(m.pos(j.vertex()));
                    Vec3d v1(m.pos(j.next().vertex()));
                    Vec3d v2(m.pos(j.opp().vertex()));
                    Vec3d v3(m.pos(j.opp().next().vertex()));
                    
                    float m1 = min(min_angle(v0, v1, v2), min_angle(v0, v2, v3));
                    float m2 = min(min_angle(v0, v1, v3), min_angle(v1, v2, v3));
                    
                    if(m1 < m2){
						// If the "cap vertex" projected onto the long edge is better in the
						// sense that there is less geometric error after the flip, then we
						// use the projected vertex. In other words, we see if it pays to snap
						// to the edge.
						Vec3d pprj = edge_dir * dot(edge_dir, p[i]-p[iplus1])+p[iplus1];
                        if(edge_error(m, he[iplus1], pprj, Vec3d(m.pos(vh[i]))) > vertex_error(m, vh[i], pprj))
                            m.pos(vh[i]) = pprj;
                        // flip if legal
                        if(precond_flip_edge(m, he[iplus1]))
                            m.flip_edge(he[iplus1]);
                        else if(boundary(m, he[iplus1]))
                            m.remove_face(*f);
                        break;
                    }
                }
            }
        }
        
    }
    
    void remove_needles(Manifold& m, float _thresh, bool averagePositions)
    {
        double thresh = _thresh * median_edge_length(m);
        for(HalfEdgeID h : m.halfedges())
            if(m.in_use(h) && length(m, h)<thresh && precond_collapse_edge(m, h))
                m.collapse_edge(h, averagePositions);
    }
    
 
    
    VertexAttributeVector<int> cluster_vertices(Manifold& m, double rad) {

        KDTree<Vec3d, VertexID> vtree;
        
        for(auto v : m.vertices())
            if(boundary(m, v))
                vtree.insert(m.pos(v), v);
        vtree.build();
        
        VertexAttributeVector<int> cluster_id(m.allocated_vertices(),-1);
        
        int cluster_ctr=0;
        for(auto v: m.vertices())
            if(boundary(m, v) && cluster_id[v] == -1)
            {
                vector<Vec3d> keys;
                vector<VertexID> vals;
                int n = vtree.in_sphere(m.pos(v), rad, keys, vals);
                
                for(int i=0;i<n;++i)
                    cluster_id[vals[i]] = cluster_ctr;
                ++cluster_ctr;
            }
        
        return cluster_id;
    }
    

    
    int stitch_mesh(Manifold& m, const VertexAttributeVector<int>& cluster_id)
    {
        map<int, vector<HalfEdgeID>> clustered_halfedges;
        for(auto v: m.vertices()) {
            HalfEdgeID h = boundary_edge(m, v);
            if(cluster_id[v] != -1 && h != InvalidHalfEdgeID)
                clustered_halfedges[cluster_id[v]].push_back(h);
        }
        int unstitched=0;
        for(auto h0 : m.halfedges())
        {
            Walker w = m.walker(h0);
            if(w.face() == InvalidFaceID)
            {
                VertexID v0 = w.opp().vertex();
                VertexID v1 = w.vertex();
                
                int cid = cluster_id[v1];
                int cid0 = cluster_id[v0];
                if(cid0 == cid) {
//                    cout << "Warning: edge endpoints in same cluster while stitching, ignoring " << endl;
                    continue;
                }
                vector<HalfEdgeID>& stitch_candidates = clustered_halfedges[cid];
                size_t i=0;
                for(;i<stitch_candidates.size(); ++i)
                {
                    HalfEdgeID h1 = stitch_candidates[i];
                    if(m.in_use(h1))
                    {
                        Walker w = m.walker(h1);
                        if(cluster_id[w.vertex()] == cluster_id[v0]) {
                            if(m.stitch_boundary_edges(h0,h1))
                                break;
                        }
                    }
                    
                }
                if(i == stitch_candidates.size())
                    ++unstitched;
            }
        }
        return unstitched;    
    }
    
    void remove_valence_two_vertices(Manifold & m)
    {
        for(VertexID v: m.vertices())
            if(valency(m,v)==2)
            {
                Walker w = m.walker(v);
                if(precond_collapse_edge(m, w.halfedge()))
                    m.collapse_edge(w.halfedge());
            }
    }
    
    int stitch_mesh(Manifold& mani, double rad)
    {
        auto cluster_id = cluster_vertices(mani, rad);
        return stitch_mesh(mani, cluster_id);
    }
    
    void stitch_more(Manifold& mani, double rad)
    {
        int unstitched, iter;
        for(iter=0;iter<2;++iter)
        {
            unstitched = stitch_mesh(mani, rad);
            vector<HalfEdgeID> hvec;
            for(HalfEdgeIDIterator h = mani.halfedges_begin(); h != mani.halfedges_end();++h)
                if(mani.walker(*h).face() == InvalidFaceID)
                    hvec.push_back(*h);
            for(size_t i=0;i<hvec.size(); ++i)
                mani.split_edge(hvec[i]);
            
        }
    }

    
    
    void close_holes(Manifold& m, int max_size)
    {
        for (auto h: m.halfedges()) {
            Walker w = m.walker(h);
            if(w.face() == InvalidFaceID)
            {
                auto h0 = w.halfedge();
                int cnt = 0;
                do {
                    ++cnt;
                    w = w.next();
                }
                while(w.halfedge() != h0 && cnt < max_size + 1);
                if(cnt <= max_size)
                    m.close_hole(h0);
            }
        }
    }
    
    void flip_orientation(Manifold& m)
    {
        vector<Vec3d> vertices;
        VertexAttributeVector<int> idvec;
        int k=0;
        for(VertexID v: m.vertices()) {
            idvec[v] = k++;
            vertices.push_back(m.pos(v));
        }
        vector<int> indices;
        vector<int> faces;
        for(FaceID f: m.faces()) {
            faces.push_back(no_edges(m, f));
			circulate_face_cw(m, f, static_cast<std::function<void(VertexID)>>([&](VertexID v){
                indices.push_back(idvec[v]);
            }));
        }
        m.clear();
        build(m, vertices.size(), vertices[0].get(), faces.size(), &faces[0], &indices[0]);
    }
    
    
    void merge_coincident_boundary_vertices(Manifold& m, double rad) {
 
        KDTree<Vec3d, VertexID> vertex_tree;
        for(auto v: m.vertices())
            if (boundary(m, v))
                vertex_tree.insert(m.pos(v), v);
        vertex_tree.build();

        for(auto v: m.vertices())
            if (boundary(m, v))
            {
                double d = rad;
                vector<Vec3d> keys;
                vector<VertexID> vals;
                int n = vertex_tree.in_sphere(m.pos(v), d, keys, vals);
                if(n>2)
                {
//                    cout << "Ambiguity: " << n << " vertices in the same spot: ";
//                    for (int i=0;i<n;++i)
//                        cout << vals[i] << " ";
//                    cout << endl;
                    continue;
                }
                if(n==2)
                {
                    VertexID v0 = vals[0];
                    VertexID v1 = vals[1];
                    m.merge_boundary_vertices(v0, v1);
                }
            }
    }
}
