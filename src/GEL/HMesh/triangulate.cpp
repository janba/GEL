/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "triangulate.h"

#include <queue>
#include <vector>
#include <iterator>
#include <cassert>

#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;

    void get_candidates(const Manifold& m, VertexID v, vector<HalfEdgeID>& candidates)
    {
        vector<VertexID> vertices;
        vector<HalfEdgeID> hedges;

        Walker wv = m.walker(v);
        for(;!wv.full_circle(); wv = wv.circulate_vertex_cw()){
            vertices.push_back(wv.vertex());
            hedges.push_back(wv.halfedge());
        }
        int N = wv.no_steps();
        vector<VertexID> vertices_check(vertices);
        assert(N == vertices.size());

        for(int i = N - 1; i >= 0; --i){
            for(Walker w = m.walker(hedges[i]).next(); w.vertex() != vertices[(i+N-1)%N]; w = w.next()){
                if(find(vertices_check.begin(), vertices_check.end(), w.vertex()) == vertices_check.end())
                    candidates.push_back(w.halfedge());       
            }
        }
    }

    float curv(const Vec3d& p, vector<Vec3d>& vec)
    {
        size_t N = vec.size();
        std::vector<Vec3d> normals;
        for(size_t i = 0; i < N; ++i){
            Vec3d norm = normalize(cross(vec[i]-p, vec[(i+1)%N]-p));
            normals.push_back(norm);
        }
        float alpha = 0;
        for(size_t i = 0; i < N; ++i)
            alpha += (vec[i]-p).length()*acos(dot(normals[i],normals[(i+1)%N]));

        return alpha;
    }

    float get_badness(const Manifold& m, VertexID v, VertexID n)
    {
        vector<HalfEdgeID> hedges;

        Walker wv = m.walker(v);
        for(; !wv.full_circle(); wv = wv.circulate_vertex_cw())
            hedges.push_back(wv.halfedge());

        vector<Vec3d> one_ring;
        vector<Vec3d> one_ring_n;
        int N = wv.no_steps();
        for(int i = N - 1; i >= 0; --i){
            for(Walker w = m.walker(hedges[i]).next(); w.vertex() != v; w = w.next()){
                one_ring.push_back(m.pos(w.vertex()));
                if(w.vertex() != n)
                    one_ring_n.push_back(m.pos(w.vertex()));
            }
        }
        return curv(m.pos(v), one_ring) - curv(m.pos(v), one_ring_n);
    }


    const CGLA::Vec3d get_normal(const Manifold& m, VertexID v)
    {

        vector<HalfEdgeID> hedges;

        Walker wv = m.walker(v);
        for(; !wv.full_circle(); wv = wv.circulate_vertex_cw())
            hedges.push_back(wv.halfedge());

        vector<Vec3d> one_ring;
        size_t N = wv.no_steps();
        for(int i = N - 1; i >= 0; --i){
            for(Walker w = m.walker(hedges[i]).next(); w.vertex() != v; w = w.next())
                one_ring.push_back(m.pos(w.vertex()));   
        }

        Vec3d norm(0);
        N = one_ring.size();
        Vec3d p0 = m.pos(v);
        for(size_t i = 0; i < N; ++i){
            Vec3d p1 = one_ring[i];
            Vec3d p2 = one_ring[(i+1) % N];
            Vec3d e0 = normalize(p1 - p0);
            Vec3d e1 = normalize(p2 - p0);
            norm += cross(e0, e1) * acos(dot(e0, e1));
        }
        return normalize(norm);
    }

    void triangulate_by_vertex_face_split(Manifold& m)
    {
        vector<FaceID> fv;
        fv.reserve(m.no_faces());
        copy(m.faces_begin(), m.faces_end(), back_inserter(fv));

        for(size_t i = 0; i < fv.size(); ++i)
            m.split_face_by_vertex(fv[i]);
    }

    void triangulate_by_edge_face_split(Manifold& m)
    {
        vector<FaceID> fv;
        fv.reserve(m.no_faces());
        copy(m.faces_begin(), m.faces_end(), back_inserter(fv));

        for(size_t i = 0; i < fv.size(); ++i)
            triangulate_face_by_edge_split(m, fv[i]);
    }

    struct PotentialEdge
    {
        int time_tag;
        float badness;
        FaceID f;
        VertexID v0, v1;
    };

    bool operator>(const PotentialEdge& e0, const PotentialEdge& e1)
    { return e0.badness>e1.badness; }

    typedef std::priority_queue<PotentialEdge,
        std::vector<PotentialEdge>,
        std::greater<PotentialEdge> > 
        PotentialEdgeQueue;

    void insert_potential_edges(const Manifold& m, VertexAttributeVector<int>& vtouched, VertexID v, PotentialEdgeQueue& pot_edges)
    {
        vector<Vec3d> one_ring;
        vector<HalfEdgeID> candidates;
        get_candidates(m, v, candidates);

        Vec3d p0 = m.pos(v);
        Vec3d norm = normal(m, v);
        int n = candidates.size();
        for(int i = 0; i < n; ++i){
            Walker w = m.walker(candidates[i]);
            VertexID v_n = w.vertex();
            Vec3d edir = normalize(m.pos(v_n) - p0);
            Vec3d norm_n = normal(m, v_n);
            float bad = sqr(dot(edir, norm));
            float bad_n = sqr(dot(edir, norm_n));

            PotentialEdge pot;

            /* So if the edge between two vertices is not orthogonal to 
            their normals, the badness is increased. Badness is also
            increased if the normals are very different. */

            pot.badness = bad+bad_n - dot(norm_n,norm);
            pot.time_tag = vtouched[v];
            pot.v0 = v;
            pot.v1 = w.vertex();
            pot.f = w.face();

            pot_edges.push(pot);
        }
    }

    void curvature_triangulate(Manifold& m)
    {
        PotentialEdgeQueue pot_edges;
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);

        // Create candidates
        for(VertexIDIterator v = m.vertices_begin(); v!= m.vertices_end(); ++v)
            insert_potential_edges(m, vtouched, *v, pot_edges);

        while(!pot_edges.empty()){
            const PotentialEdge& pot_edge = pot_edges.top();
            // Record all the vertices of the face. We need to 
            // recompute the candidates.
            std::vector<VertexID> reeval_vec;

            for(Walker w = m.walker(pot_edge.f); !w.full_circle(); w = w.circulate_face_cw())
                reeval_vec.push_back(w.vertex());

            // Make sure that the vertex has not been touched since 
            // we created the potential edge. If the vertex has been
            // touched the candidate edge has either (a) been processed,
            // (b) received a lower priority or (c) become invalid.
            if(pot_edge.time_tag == vtouched[pot_edge.v0]){
                vector<Vec3d> one_ring;
                vector<HalfEdgeID> candidates;

                m.split_face_by_edge(pot_edge.f, pot_edge.v0, pot_edge.v1);

                // Recompute priorities.
                for(size_t i = 0; i < reeval_vec.size(); ++i){
                    VertexID& v = reeval_vec[i];
                    ++vtouched[v];
                    insert_potential_edges(m, vtouched, v, pot_edges);
                }

            }
            pot_edges.pop();
        }

    }

    void shortest_edge_triangulate(Manifold& m)
    {
        int work;
        do{
            // For every face.
            work = 0;
            for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
                // Create a vector of vertices.
                vector<VertexID> verts;
                for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
				{
					FaceID fa = w.face();
					FaceID fb = *f;
					assert(fa==fb);
                    verts.push_back(w.vertex());
				}
                // If there are just three we are done.
                if(verts.size() == 3) continue;

                // Find vertex pairs that may be connected.
                vector<pair<int,int> > vpairs;
                const int N = verts.size();
                for(int i = 0; i < N - 2; ++i){
                    for(int j = i + 2; j < N; ++j){
                        if(verts[i] != verts[j] && !connected(m, verts[i], verts[j]))
                            vpairs.push_back(pair<int,int>(i, j));
                    }
                }
                if(vpairs.empty()){
                    cout << "Warning: could not triangulate a face." 
                        << "Probably a vertex appears more than one time in other vertex's one-ring" << endl;
                    continue;
                }

                /* For all vertex pairs, find the edge lengths. Combine the
                vertices forming the shortest edge. */

                float min_len=FLT_MAX;
                int min_k = -1;
                for(size_t k = 0; k < vpairs.size(); ++k){
                    int i = vpairs[k].first;
                    int j = vpairs[k].second;
                    float len = sqr_length(m.pos(verts[i]) - m.pos(verts[j]));

                    if(len<min_len){
                        min_len = len;
                        min_k = k;
                    }
                }
                assert(min_k != -1);

                {
                    // Split faces along edge whose midpoint is closest to isovalue
                    int i = vpairs[min_k].first;
                    int j = vpairs[min_k].second;
                    m.split_face_by_edge(*f, verts[i], verts[j]);
                }
                ++work;

            }
        }
        while(work);
    }

    void triangulate_face_by_edge_split(Manifold& m, FaceID f)
    {
        if(no_edges(m, f)<=3)
            return;
        
        Walker w = m.walker(f);

        // as long as f is not a triangle
        while(w.next().next().next().halfedge() != w.halfedge()){
            // assert that face has at least 3 edges
            // f is split into triangle from first three vertices, and becomes remaining polygon in next iteration
            assert(w.next().next().halfedge() != w.halfedge());
            VertexID v0 = w.vertex();
            VertexID v1 = w.next().next().vertex();
            FaceID f_old = f;
            if(v0 != v1 && !connected(m, v0, v1))
                f = m.split_face_by_edge(f, v0, v1);
            if(f == f_old)
                return;
        }
    }
}
