/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/dual.h>
#include <GEL/HMesh/subdivision.h>

#include <vector>
#include <GEL/CGLA/Vec3d.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    
    void loop_split(Manifold& m_in, Manifold& m)
    {
        if(&m != &m_in)
            m = m_in;
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);
        
        vector<HalfEdgeID> hedges;
        for(HalfEdgeID h : m.halfedges())
            if(h<m.walker(h).opp().halfedge())
                hedges.push_back(h);
        
        for(HalfEdgeID h : hedges)
            vtouched[m.split_edge(h)] = 1;
        
        FaceSet faces(m.faces());
        for(FaceID fid : faces)
        {
            Walker w = m.walker(fid);
            
            if(vtouched[w.vertex()] == 0)
                w = w.next();
            
            assert(vtouched[w.vertex()] == 1);
            
            VertexID v1, orig_vert = w.vertex();
            w = w.next();
            FaceID f = fid;
            do
            {
                VertexID v0 = w.opp().vertex();
                w = w.next();
                v1  = w.vertex();
                w = w.next();
                assert(vtouched[v0] == 1);
                assert(vtouched[v1] == 1);
                f = m.split_face_by_edge(f, v0, v1);
            }
            while (v1 != orig_vert);
        }
    }
    
    void cc_split(Manifold& m, Manifold& m_out)
    {
        VertexAttributeVector<int> touched(m.no_vertices(), -1);
        // Mark all original vertices
        for (auto v: m.vertices())
            touched[v] = 1;
        
        // Split all of the original edges by inserting
        // a vertex on them.
        vector<HalfEdgeID> orig_edges;
        for (auto h: m.halfedges())
            if (h == m.walker(h).hmin())
                orig_edges.push_back(h);
        for (auto h: orig_edges)
            m.split_edge(h);
        
        // Split all vertices by inserting a vertex in the middle
        vector fvec(begin(m.faces()), end(m.faces()));
        for(auto f: fvec) {
            auto vnew = m.split_face_by_vertex(f);
            touched[vnew] = 2;
        }

        // Remove edges from the face points to original vertices.
        for(auto h: m.halfedges())
            if (m.in_use(h)) {
                auto w = m.walker(h);
                if(touched[w.vertex()] == 1 && touched[w.opp().vertex()] == 2)
                    m.merge_faces(w.face(), h);
            }
        m_out = m;
    }
    
    void root3_subdivide(Manifold& m_in, Manifold& m)
    {
        if(&m != &m_in)
            m = m_in;
        
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);
        VertexAttributeVector<Vec3d> new_pos(m.allocated_vertices(), Vec3d(0));
        
        for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid) {
            int v = valency(m, *vid);
            double beta = (4.0 - 2.0 * cos(2.0 * M_PI / v))/(9.0*v);
            new_pos[*vid] = (1.0 - v * beta) * m.pos(*vid);
            for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                new_pos[*vid] += beta * m.pos(w.vertex());
            }
        }

        vector<FaceID> faces;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
            faces.push_back(*f);
        for(int i=0;i<faces.size(); ++i)
            vtouched[m.split_face_by_vertex(faces[i])] = 1;
        
        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h)
        {
            Walker w = m.walker(*h);
            
            if(vtouched[w.vertex()] == 0 && vtouched[w.opp().vertex()] == 0 &&
                precond_flip_edge(m, *h))
                m.flip_edge(*h);
        }
        
        for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            if(vtouched[*vid] == 0)
                m.pos(*vid) = new_pos[*vid];
    }

    void rootCC_subdivide(Manifold& m_in, Manifold& m)
    {
        if(&m != &m_in)
            m = m_in;
        
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);
        VertexAttributeVector<Vec3d> new_pos(m.allocated_vertices(), Vec3d(0));
        
        for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid) {
            int v = valency(m, *vid);
            double beta = (4.0 - 2.0 * cos(2.0 * M_PI / v))/(9.0*v);
            new_pos[*vid] = (1.0 - v * beta) * m.pos(*vid);
            for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                new_pos[*vid] += beta * m.pos(w.vertex());
            }
        }
        
        vector<FaceID> faces;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
            faces.push_back(*f);
        for(int i=0;i<faces.size(); ++i)
            vtouched[m.split_face_by_vertex(faces[i])] = 1;
        
        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h)
        {
            Walker w = m.walker(*h);
            
            if(vtouched[w.vertex()] == 0 && vtouched[w.opp().vertex()] == 0)
                m.merge_faces(w.face(), *h);
        }
        
        for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            if(vtouched[*vid] == 0)
                m.pos(*vid) = new_pos[*vid];
    }

    
    void butterfly_subdivide(Manifold& m_in, Manifold& m)
    {
        const float S[4][6] = {{5.0/12.0, -1.0/12.0, -1.0/12.0, 0, 0, 0},
            {3.0/8.0, 0, -1.0/8.0, 0, 0, 0},
            {0.35,
                0.03090169943749475,
                -0.08090169943749474,
                -0.08090169943749474,
                0.03090169943749468,0 },
            {0, 1.0f/8, -1.0f/8, 0, -1.0f/8, 1.0f/8}};

        
        if(&m != &m_in)
            m = m_in;

        HalfEdgeAttributeVector<Vec3d> new_vertices_pos(m.allocated_halfedges(), Vec3d(0.0));
        HalfEdgeAttributeVector<int> htouched(m.allocated_halfedges(), 0);
        VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);

        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
            {
                Walker w = m.walker(*hid);
                VertexID v0 = w.opp().vertex();
                
                int K = valency(m, v0);
                int Ko = valency(m, w.vertex());
                if((K==6 && Ko==6)|| K != 6)
                {
                    Vec3d pos((K==6 ? 1.0 : 3.0/4.0) * m.pos(v0));
                    for(int k=0;k<K; ++k, w = w.circulate_vertex_ccw())
                    {
                        double s = (K<=6) ? S[K-3][k]:(0.25+cos((2.0*M_PI*k)/K)+0.5*cos((4.0*M_PI*k)/K))/K;
                        pos += s * m.pos(w.vertex());                        
                    }
                    new_vertices_pos[*hid] = pos;
                    htouched[*hid] = 1;
                }
            }
        loop_split(m, m);

        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
            if(htouched[*hid])
            {
                VertexID v = m.walker(*hid).opp().vertex();
                if(vtouched[v] == 0)
                    m.pos(v) =Vec3d(0);
                m.pos(v) += new_vertices_pos[*hid];
                vtouched[v]+=1;
            }
        for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
            if(vtouched[*vid])
                m.pos(*vid) /= vtouched[*vid];
    }
    
    enum Subd {QUAD_SUBD, CC_SUBD, LOOP_SUBD, TRI_SUBD};
    
    void subd_smooth(Subd subd_method, Manifold& m)
    {
        VertexAttributeVector<Vec3d> new_vertices(Vec3d(0));
        
        for(auto f: m.faces())
        {
            circulate_face_ccw(m, f, [&](VertexID v0) {
                double val = valency(m, v0);
                double A,B;
                
                switch(subd_method)
                {
                    case QUAD_SUBD:
                        A = 1.0 / (4.0 * val);
                        B = 1.0 / (4.0 * val);
                        break;
                    case CC_SUBD:
                        A = (1.0-3.0/val) * (1.0/val);
                        B = sqr(1.0/val);
                        break;
                    case TRI_SUBD:
                        A = 2.0 / (8.0 * val);
                        B = 3.0 / (8.0 * val);
                        break;
                    case LOOP_SUBD:
                        float w = 5.0/8.0 - sqr(3.0/8.0 + 0.25 * cos(2.0*M_PI/val));
                        A = (1.0-2.0*w)/val;
                        B = w/val;
                        break;
                }
                circulate_face_ccw(m, f, [&](VertexID v) {
                    if(v == v0)
                        new_vertices[v0] += A * m.pos(v);
                    else
                        new_vertices[v0] += B * m.pos(v);
                });
                
            });
        }
        m.positions = new_vertices;
    }

    void cc_smooth(Manifold& m)
    {
        subd_smooth(CC_SUBD, m);
    }

    void volume_preserving_cc_smooth(Manifold& m, double w, int iter) {
        FaceAttributeVector<double> face_area;
        for (auto f: m.faces())
            face_area[f] = area(m, f);
        VertexAttributeVector<double> one_ring_area(m.allocated_vertices(),0);
        double avg_one_ring_area = 0;
        for (auto v: m.vertices()) {
            for (auto f: m.incident_faces(v))
                one_ring_area[v] += face_area[f];
            avg_one_ring_area += one_ring_area[v];
        }
        avg_one_ring_area /= m.no_vertices();
        for (int i=0;i<iter;++i) {
            auto non_smooth_pos = m.positions_attribute_vector();
            subd_smooth(CC_SUBD, m);
            for (auto v: m.vertices()) {
                double w = min(0.95, 0.25*one_ring_area[v]/avg_one_ring_area);
                m.pos(v) = (1-w) * m.pos(v) + w * non_smooth_pos[v];
            }
        }

//        for (int i=0;i<iter;++i) {
//            subd_smooth(CC_SUBD, m);
//            auto old_pos = m.positions_attribute_vector();
//            subd_smooth(CC_SUBD, m);
//            for (auto v: m.vertices()) {
//                double w_local = w*(one_ring_area[v]/max_one_ring_area);
//                m.pos(v) = old_pos[v] - w_local *(m.pos(v) - old_pos[v]);
//            }
//        }
    }

    
    void loop_smooth(Manifold& m)
    {
        subd_smooth(LOOP_SUBD, m);
    }


}
