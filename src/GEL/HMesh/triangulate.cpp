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

    FaceID clip_ear(HMesh::Manifold& m, FaceID f) {
        int N = no_edges(m, f);
        double total_area = area(m,f);
        double ideal_ear_area = total_area / (N-2);
        if(N<4)
            return InvalidFaceID;
        
        Vec3d norm = normal(m,f);

        double max_energy = -1.0;
        HalfEdgeID h_max = InvalidHalfEdgeID;
        circulate_face_ccw(m,f,[&](Walker w){
            VertexID v0 = w.opp().vertex();
            VertexID v1 = w.next().vertex();
            if(v0 != v1 && !connected(m, v0, v1)) {
                Vec3d p = m.pos(w.vertex());
                Vec3d pp = (m.pos(v0) - p);
                Vec3d pn = (m.pos(v1) - p);
                Vec3d area_vec = cross(pn,pp);
                double ear_area = 0.5 * length(area_vec);
                double area_energy = 1.0-min(1.0,max(0.0,(ear_area-ideal_ear_area)/(total_area-ideal_ear_area)));
                double convexity = dot(norm, area_vec);
                double energy = area_energy*dot(pn,pp)/(length(pp)*length(pn));
                if (convexity>0.0 && energy>max_energy)
                {
                    max_energy = energy;
                    h_max=w.halfedge();
                }
            }
        });
        if(h_max != InvalidHalfEdgeID)
        {
            Walker w = m.walker(h_max);
            FaceID f_new = m.split_face_by_edge(f, w.opp().vertex(), w.next().vertex());
            if(f_new != InvalidFaceID)
                return f_new;
        }
        return InvalidFaceID;
    }

    FaceID shortest_edge_split(Manifold& m, FaceID f) {
        // Create a vector of vertices.
        vector<VertexID> verts;
        int no_edges = circulate_face_ccw(m, f, [&verts](VertexID vn){
            verts.push_back(vn);
        });
        // If there are just three we are done.
        if(no_edges == 3) return InvalidFaceID;
        
        using VP = pair<int,int>;
        // Find vertex pairs that may be connected.
        vector<VP> vpairs;
        const int N = verts.size();
        for(int i = 0; i < N - 2; ++i){
            for(int j = i + 2; j < N; ++j){
                if(verts[i] != verts[j] && !connected(m, verts[i], verts[j]))
                    vpairs.push_back(pair<int,int>(i, j));
            }
        }
        if(! vpairs.empty()){
            /* For all vertex pairs, find the edge lengths. Combine the
             vertices forming the shortest edge. */
            auto mei = min_element(begin(vpairs), end(vpairs), [&](const VP& a, const  VP& b) {
                double lena = sqr_length(m.pos(verts[a.first]) - m.pos(verts[a.second]));
                double lenb = sqr_length(m.pos(verts[b.first]) - m.pos(verts[b.second]));
                return lena < lenb;
            });
            int i = mei->first;
            int j = mei->second;
            FaceID f_new = m.split_face_by_edge(f, verts[i], verts[j]);
            return f_new;
        }
        return InvalidFaceID;
    }

    int triangulate(Manifold& m, FaceID f, TriangulationMethod policy) {
        queue<FaceID> fq;
        fq.push(f);
        int work = 0;
        while(!fq.empty()) {
            FaceID f = fq.front();
            fq.pop();
            FaceID fn = policy==CLIP_EAR ? clip_ear(m, f) : shortest_edge_split(m, f);
            if (fn != InvalidFaceID) {
                fq.push(fn);
                if(policy == SHORTEST_EDGE)
                    fq.push(f);
                ++work;
            }
        }
        return work;
    }

    void triangulate(Manifold& m, TriangulationMethod policy)
    {
        int work;
        do{
            // For every face.
            work = 0;
            for (auto f: m.faces()) {
                work += triangulate(m, f, policy);
            }
        }
        while(work);
    }

}
