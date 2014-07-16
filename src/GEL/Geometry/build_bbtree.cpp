/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "build_bbtree.h"

#include "../HMesh/Manifold.h"

using namespace CGLA;
using namespace std;
using namespace HMesh;

namespace
{
    const float EDGE_MIN_SQ_LENGTH = CGLA::MINUTE;

    inline bool degenerate_edge(const Manifold& m, HalfEdgeID h)
    {
		Walker w = m.walker(h);
        if(sqr_length(m.pos(w.vertex()) - m.pos(w.opp().vertex())) < 1e-8)
            return true;
        return false;
    }

}

namespace Geometry
{

    template<class BBTree>
    void build_tree_robust(Manifold& m, BBTree& tree)
    {
        vector<Triangle> triangle_vec;

        for(FaceIDIterator fi=m.faces_begin(); fi != m.faces_end();++fi)
        {
            Vec3d face_normal = normal(m, *fi);
            
            Walker w = m.walker(*fi);

            Vec3f v0,v1,v2;
            Vec3f vn0,vn1,vn2;
            Vec3f en0,en1,en2;

            v0  = Vec3f(m.pos(w.vertex()));
            vn0 = Vec3f(normal(m, w.vertex()));
            FaceID adj_f =  w.next().opp().face();
            if(adj_f != InvalidFaceID)
                en0 = Vec3f(normalize(face_normal + normal(m, adj_f)));
            else
                en0 = Vec3f(face_normal);
            
            w = w.next();

            v1  = Vec3f(m.pos(w.vertex()));
            vn1 = Vec3f(normal(m, w.vertex()));
            adj_f =  w.next().opp().face();
            if(adj_f != InvalidFaceID)
                en1 = Vec3f(normalize(face_normal + normal(m, adj_f)));
            else
                en1 = Vec3f(face_normal);

            
            w = w.next();

            v2  = Vec3f(m.pos(w.vertex()));
            vn2 = Vec3f(normal(m, w.vertex()));
            adj_f =  w.next().opp().face();
            if(adj_f != InvalidFaceID)
                en2 = Vec3f(normalize(face_normal + normal(m, adj_f)));
            else
                en2 = Vec3f(face_normal);
            
            if(sqr_length(v0-v1)>EDGE_MIN_SQ_LENGTH &&
                sqr_length(v1-v2)>EDGE_MIN_SQ_LENGTH &&
                sqr_length(v2-v0)>EDGE_MIN_SQ_LENGTH)
                triangle_vec.push_back(Triangle(v0,v1,v2,vn0,vn1,vn2,en0,en1,en2));
            else
                cout << "Killing degenerate triangle" << endl;
        }

        tree.build(triangle_vec);
    }

    void build_OBBTree(HMesh::Manifold& m, OBBTree& tree)
    {
        build_tree_robust<OBBTree>(m, tree);
    }

    void build_AABBTree(HMesh::Manifold& m, AABBTree& tree)
    {
        build_tree_robust<AABBTree>(m, tree);
    }

}
