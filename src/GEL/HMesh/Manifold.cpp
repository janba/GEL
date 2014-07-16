/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */
#include <iterator>
#include "Manifold.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <iterator>

#include "../Geometry/TriMesh.h"

namespace HMesh
{
	
    using namespace std;
    using namespace Geometry;
    using namespace CGLA;
	
    namespace
    {
        /************************************************************
		 * Edgekeys and Edges for halfedge identification during build
		 ************************************************************/
        class EdgeKey
        {
        public:
            EdgeKey(VertexID va, VertexID vb)
            {
                if(va < vb){
                    v0 = va;
                    v1 = vb;
                }
                else{
                    v0 = vb;
                    v1 = va;
                }
            }
			
            bool operator<(const EdgeKey& k2) const
            {
                if(v0 < k2.v0){
                    return true;
                }
                else if( k2.v0 < v0){
                    return false;
                }
                else{
                    return v1 < k2.v1;
                }
            }
            
            bool operator==(const EdgeKey& k2) const {return k2.v0 == v0 && k2.v1==v1;}
            
            size_t hash() const {return v0.get_index()*3125*49+v1.get_index()*3125+7;}
        private:
            VertexID v0;
            VertexID v1;
        };
		
        struct Edge
        {
            HalfEdgeID h0;
            HalfEdgeID h1;
            int count;
            Edge() : count(0){}
        };
    }
	
    /*********************************************
	 * Public functions
	 *********************************************/
    void Manifold::build(const TriMesh& mesh)
    {
        // A vector of 3's - used to tell build how many indices each face consists of
        vector<int> faces(mesh.geometry.no_faces(), 3);
		
        build_template( static_cast<size_t>(mesh.geometry.no_vertices()),
					   reinterpret_cast<const float*>(&mesh.geometry.vertex(0)),
					   static_cast<size_t>(faces.size()),
					   static_cast<const int*>(&faces[0]),
					   reinterpret_cast<const int*>(&mesh.geometry.face(0)));
    }
    
    void Manifold::build(   size_t no_vertices,
						 const float* vertvec,
						 size_t no_faces,
						 const int* facevec,
						 const int* indices)
    {
        build_template(no_vertices, vertvec, no_faces, facevec, indices);
    }
    
    void Manifold::build(   size_t no_vertices,
						 const double* vertvec,
						 size_t no_faces,
						 const int* facevec,
						 const int* indices)
    {
        build_template(no_vertices, vertvec, no_faces, facevec, indices);
    }
    
    FaceID Manifold::add_face(std::vector<Manifold::Vec> points)
    {
        int F = points.size();
        vector<int> indices;
        for(size_t i=0;i<points.size(); ++i)
            indices.push_back(i);
        FaceID fid = *faces_end();
        build(points.size(), reinterpret_cast<double*>(&points[0]), 1, &F, &indices[0]);
        return fid;
    }
    
    bool Manifold::remove_face(FaceID fid)
    {
        if(!in_use(fid))
            return false;
        
        HalfEdgeID he = kernel.last(fid);
        HalfEdgeID last = he;
        
        vector<HalfEdgeID> halfedges;
        do
        {
            halfedges.push_back(he);
            kernel.set_face(he, InvalidFaceID);
            he = kernel.next(he);
        }
        while(he != last);
        
        vector<HalfEdgeID> halfedges_garbage;
        vector<VertexID> vertices_garbage;
        for(size_t i=0;i<halfedges.size(); ++i)
        {
            HalfEdgeID h = halfedges[i];
            HalfEdgeID hopp = kernel.opp(h);
            if(kernel.face(hopp) == InvalidFaceID)
            {
                halfedges_garbage.push_back(h);
                VertexID v0 = kernel.vert(hopp);
                if(valency(*this, v0) <= 2 )
                    vertices_garbage.push_back(v0);
                else {
                    link(kernel.prev(h), kernel.next(hopp));
                    kernel.set_out(v0, kernel.opp(kernel.prev(h)));
                }
                VertexID v1 = kernel.vert(h);
                if(valency(*this, v1)>2)
                {
                    link(kernel.prev(hopp),kernel.next(h));
                    kernel.set_out(v1, kernel.opp(kernel.prev(hopp)));
                }
                
            }
        }
        for(size_t i=0;i<halfedges_garbage.size();++i){
            kernel.remove_halfedge(kernel.opp(halfedges_garbage[i]));
            kernel.remove_halfedge(halfedges_garbage[i]);
        }
        for(size_t i=0;i<vertices_garbage.size(); ++i)
            kernel.remove_vertex(vertices_garbage[i]);
        
        kernel.remove_face(fid);
        return true;
    }
    

    bool Manifold::remove_edge(HalfEdgeID hid)
    {
        if(!in_use(hid))
            return false;
        
        FaceID f0 = kernel.face(hid);
        FaceID f1 = kernel.face(kernel.opp(hid));
        
        remove_face(f0);
        remove_face(f1);

        return true;
    }

    
    
    bool Manifold::remove_vertex(VertexID vid)
    {
        if(!in_use(vid))
            return false;
    
        vector<FaceID> faces;
        int N = circulate_vertex_ccw(*this, vid, (std::function<void(FaceID)>)[&](FaceID f) {
            faces.push_back(f);
        });
        for(size_t i=0;i<N;++i)
            remove_face(faces[i]);
            
        return true;
    }


	
    void Manifold::collapse_edge(HalfEdgeID h, bool avg_vertices)
    {
        HalfEdgeID ho = kernel.opp(h);
        VertexID hv = kernel.vert(h);
        VertexID hov = kernel.vert(ho);
        HalfEdgeID hn = kernel.next(h);
        HalfEdgeID hp = kernel.prev(h);
        HalfEdgeID hon = kernel.next(ho);
        HalfEdgeID hop = kernel.prev(ho);
		FaceID f = kernel.face(h);
		FaceID fo = kernel.face(ho);
		
        // average the vertex positions
        pos(hv) = avg_vertices ? (0.5f * (pos(hov) + pos(hv))) : pos(hv);
		
        // update all halfedges pointing to hov to point to hv, effectively removing hov from all loops
        HalfEdgeID he = kernel.out(hov);
        HalfEdgeID last = he;
		do {
			assert(kernel.vert(kernel.opp(he)) == hov);
            kernel.set_vert(kernel.opp(he), hv);
            he = kernel.next(kernel.opp(he));
        }
		while(he != last);
        kernel.set_out(hv, hn);
		
        // link hp and hn, effectively removing h from opposite loop
        link(hp, hn);
        // make face owning h own hn instead
        if(kernel.face(h) != InvalidFaceID)
            kernel.set_last(f, hn);
		
        // link hop and hon, effectively removing h from opposite face loop
        link(hop, hon);
        // make opposite face owning h own hon instead
        if(kernel.face(ho) != InvalidFaceID)
            kernel.set_last(fo, hon);
		
        // remove the obsolete entities
        kernel.remove_vertex(hov);
        kernel.remove_halfedge(h);
        kernel.remove_halfedge(ho);
		
        // verify that remaining faces haven't become degenerate because of collapse
        remove_face_if_degenerate(hn);
        remove_face_if_degenerate(hon);
		
        // verify that v is sane after collapse
        ensure_boundary_consistency(hv);
    }
	
    FaceID Manifold::split_face_by_edge(FaceID f, VertexID v0, VertexID v1)
    {
		if(connected(*this, v0, v1))
            return InvalidFaceID;
        
        HalfEdgeID he = kernel.last(f);
        HalfEdgeID last = he;
        int steps = 0;
        while(he != last || steps == 0){
            ++steps;
            he = kernel.next(he);
        }
		
        // make sure this is not a triangle
        assert(steps > 3);
        // make sure we are not trying to connect a vertex to itself
        assert(v0 != v1);
		
		HalfEdgeID h0 = kernel.out(v0);
		for(Walker w = walker(v0); !w.full_circle(); w = w.circulate_vertex_cw()){
			if(w.face() == f){
				h0 = w.halfedge();
				break;
			}
		}
		assert(kernel.face(h0) != InvalidFaceID);
        assert(kernel.face(h0) == f);
		
        // the halfedge belonging to f, going out from v0, is denoted h. Move along h until we hit v1.
        // now we have the halfedge which belongs to f and points to v1.
        HalfEdgeID h = h0;
        while(kernel.vert(h) != v1){
            h = kernel.next(h);
        }
        assert(h != h0);
		
        // create a new halfedge ha which connects v1 and v0 closing the first loop.
        HalfEdgeID h1 = kernel.next(h);
        HalfEdgeID ha = kernel.add_halfedge();
        link(h, ha);
        link(ha, h0);
        kernel.set_face(ha, f);
        kernel.set_vert(ha, v0);
        kernel.set_last(f, ha);
		
        // create a new face, f2, and set all halfedges in the remaining part of the polygon to point to this face.
        h = h1;
        FaceID f2 = kernel.add_face();
        while(kernel.vert(h) != v0){
            kernel.set_face(h, f2);
            h = kernel.next(h);
        }
        kernel.set_face(h, f2);
        assert(h != h1);
		
        // create a new halfedge hb to connect v0 and v1.
        HalfEdgeID hb = kernel.add_halfedge();
        link(h, hb);
        link(hb, h1);
        kernel.set_face(hb, f2);
        kernel.set_vert(hb, v1);
        kernel.set_last(f2, hb);
		
        // complete the operation by gluing the two new halfedges
        glue(ha, hb);
		
        // assert sanity of operation
        assert(kernel.next(kernel.opp(kernel.prev(h1))) == h0);
        assert(kernel.next(kernel.opp(kernel.prev(h0))) == h1);
        assert(kernel.face(kernel.next(hb)) == f2);
        assert(kernel.face(kernel.next(kernel.next(hb))) == f2);
        assert(kernel.face(hb) == f2);
		
        // return handle to newly created face
        return f2;
    }
	
    VertexID Manifold::split_face_by_vertex(FaceID f)
    {
        //create the new vertex, with the barycenter of the face as position
        Manifold::Vec p(0.0f);
        HalfEdgeID last_he = kernel.last(f);
        HalfEdgeID he = last_he;
        int steps = 0;
		
        while(he != last_he || steps == 0){
            p += positions[kernel.vert(he)];
            ++steps;
            he = kernel.next(he);
        }
		
        p /= steps;
		
        VertexID v = kernel.add_vertex();
        positions[v]  = p;
		
        //circulate the face, create halfedges and connect to vertex
        vector<HalfEdgeID> eout;
        vector<HalfEdgeID> ein;
        last_he = kernel.last(f);
        he = last_he;
		
        do{
            HalfEdgeID hn = kernel.next(he);
			
            HalfEdgeID ho = kernel.add_halfedge();
            HalfEdgeID hi = kernel.add_halfedge();
			
            FaceID fn = kernel.add_face();
            kernel.set_face(hi, fn);
            kernel.set_vert(hi, v);
            kernel.set_face(ho, fn);
            kernel.set_vert(ho, kernel.vert(kernel.opp(he)));
            kernel.set_face(he, fn);
            kernel.set_last(fn, ho);
			
            link(hi, ho);
            link(ho, he);
            link(he, hi);
			
            eout.push_back(ho);
            ein.push_back(hi);
			
            he = hn;
        }
        while(he != last_he);
		
        // glue new halfedges together
        size_t N = eout.size();
        for(size_t i = 0; i < N; ++i){
            glue(ein[i], eout[(i+1)%N]);
        }
        kernel.set_out(v, eout[0]);
		
        //remove the now replaced face
        kernel.remove_face(f);
		
        return v;
    }
    VertexID Manifold::split_edge(HalfEdgeID h)
    {
        HalfEdgeID ho = kernel.opp(h);
        VertexID v = kernel.vert(h);
        VertexID vo = kernel.vert(ho);
		
        //create the new vertex with middle of edge as position and update connectivity
        VertexID vn = kernel.add_vertex();
        positions[vn] = .5f * (positions[v] + positions[vo]);
        kernel.set_out(vn, h);
		
        //create two necessary halfedges, and update connectivity
        HalfEdgeID hn = kernel.add_halfedge();
        HalfEdgeID hno = kernel.add_halfedge();
		
        kernel.set_out(vo, hn);
        kernel.set_out(v, hno);
		
        glue(h, hno);
        glue(hn, ho);
		
        link(kernel.prev(h), hn);
        link(hn, h);
        kernel.set_vert(hn, vn);
		
        link(kernel.prev(ho), hno);
        link(hno, ho);
        kernel.set_vert(hno, vn);
		
        // update faces in case of non boundary edge
        if(kernel.face(h) != InvalidFaceID)
            kernel.set_last(kernel.face(h), hn);
		
        if(kernel.face(ho) != InvalidFaceID)
            kernel.set_last(kernel.face(ho), ho);
		
        // update new edge with faces from existing edge
        kernel.set_face(hn, kernel.face(h));
        kernel.set_face(hno, kernel.face(ho));
		
        //if split occurs on a boundary, consistency must be ensured.
        ensure_boundary_consistency(vn);
        ensure_boundary_consistency(v);
        ensure_boundary_consistency(vo);
		
        return vn;
    }
    
    size_t link_intersection(const Manifold& m, VertexID v0, VertexID v1, vector<VertexID>& lisect)
    {
        // get the one-ring of v0
        vector<VertexID> link0;
        circulate_vertex_ccw(m, v0, (std::function<void(VertexID)>)[&](VertexID vn) {
            link0.push_back(vn);
        });
		
        // get the one-ring of v1
        vector<VertexID> link1;
        circulate_vertex_ccw(m, v1, (std::function<void(VertexID)>)[&](VertexID vn) {
            link1.push_back(vn);
        });
		
        // sort the vertices of the two rings
        sort(link0.begin(), link0.end());
        sort(link1.begin(), link1.end());
		
        // get the intersection of the shared vertices from both rings
        std::back_insert_iterator<vector<VertexID> > lii(lisect);
        set_intersection(link0.begin(), link0.end(),
						 link1.begin(), link1.end(),
						 lii);
        
        return lisect.size();
    }
    
    bool Manifold::stitch_boundary_edges(HalfEdgeID h0, HalfEdgeID h1)
    {
        // Cannot stitch an edge with itself
        if(h0 == h1)
            return false;
        
        // Only stitch a pair of boundary edges.
        if(kernel.face(h0) == InvalidFaceID && kernel.face(h1) == InvalidFaceID)
        {
            HalfEdgeID h0o = kernel.opp(h0);
            HalfEdgeID h1o = kernel.opp(h1);
            VertexID v0a = kernel.vert(h0);
            VertexID v0b = kernel.vert(kernel.opp(h1));
            VertexID v1b = kernel.vert(h1);
            VertexID v1a = kernel.vert(kernel.opp(h0));
            
            // Case below implies that h0 and h1 are the same edge with different ID
            // That should not happen.
            assert(!(v0a == v0b && v1a == v1b));
            
            if(connected(*this, v0a, v0b))
                return false;
            if(connected(*this, v1a, v1b))
                return false;
            
            
            if(v0b != v0a)
            {
                
                // Check the link intersection v0a and v0b are welded together
                // if they share a neighbouring vertex, it will appear twice in the combined
                // one ring unless it v1a and v1a==v1b
                vector<VertexID> lisect;
                if(link_intersection(*this, v0a, v0b, lisect))
                {
                    vector<VertexID>::iterator iter;
                    
                    if(v1a == v1b)
                    {
                        iter = find(lisect.begin(), lisect.end(), v1a);
                        if(iter != lisect.end())
                            lisect.erase(iter);
                    }
                    iter = find(lisect.begin(), lisect.end(), kernel.vert(kernel.next(h0)));
                    if(iter != lisect.end())
                        lisect.erase(iter);
                    if(lisect.size() > 0)
                        return false;
                }
            }
            
            if(v1a != v1b)
            {
                // Check the same for the other endpoints.
                vector<VertexID> lisect;
                if(link_intersection(*this, v1a, v1b, lisect))
                {
                    vector<VertexID>::iterator iter;
                    
                    if(v0a == v0b)
                    {
                        iter = find(lisect.begin(), lisect.end(), v0a);
                        if(iter != lisect.end())
                            lisect.erase(iter);
                    }
                    iter = find(lisect.begin(), lisect.end(), kernel.vert(kernel.next(h1)));
                    if(iter != lisect.end())
                        lisect.erase(iter);
                    if(lisect.size() > 0)
                        return false;
                }
            }
            
            
            if(v0b != v0a)
                circulate_vertex_ccw(*this, v0b, (std::function<void(Walker&)>)[&](Walker hew) {
                    kernel.set_vert(hew.opp().halfedge(), v0a);
                });
            
            if(v1b != v1a)
                circulate_vertex_ccw(*this, v1b, (std::function<void(Walker&)>)[&](Walker hew) {
                    kernel.set_vert(hew.opp().halfedge(), v1a);
                });
            
            if(v0a != v0b)
            {
                HalfEdgeID h1p = kernel.prev(h1);
                HalfEdgeID h0n = kernel.next(h0);
                
                if(kernel.next(h0n) == h1p)
                {
                    glue(kernel.opp(h0n), kernel.opp(h1p));
                    kernel.set_out(kernel.vert(h0n),kernel.opp(h0n));
                    kernel.remove_halfedge(h0n);
                    kernel.remove_halfedge(h1p);
                }
                else
                    link(h1p, h0n);
                kernel.remove_vertex(v0b);
            }
            
            if(v1a != v1b)
            {
                HalfEdgeID h0p = kernel.prev(h0);
                HalfEdgeID h1n = kernel.next(h1);
                if(kernel.next(h1n) == h0p)
                {
                    glue(kernel.opp(h0p), kernel.opp(h1n));
                    kernel.set_out(kernel.vert(h1n), kernel.opp(h1n));
                    kernel.remove_halfedge(h0p);
                    kernel.remove_halfedge(h1n);
                }
                else
                    link(h0p, h1n);
                kernel.remove_vertex(v1b);
            }
            glue(h0o, h1o);
            
            kernel.remove_halfedge(h0);
            kernel.remove_halfedge(h1);
            
            kernel.set_out(v1a, h1o);
            kernel.set_out(v0a, h0o);
            
            ensure_boundary_consistency(v1a);
            ensure_boundary_consistency(v0a);
            
            return true;
        }
        return false;
    }
    
    
    
    FaceID Manifold::merge_one_ring(VertexID v, float max_loop_length)
    {
        // If the vertex is either not in use or has just
        // one incident edge (or less), bail out.
        int val = valency(*this,v);
        if(!in_use(v) || val<2)
            return InvalidFaceID;
        
        // If the vertex is  a boundary vertex, we preparte the walker so that the
        // first face visited is not the invalid face outside the boundary. If the boundary
        // vertex is adjacent to only one vertex, there is little to do and we bail out.
        bool vertex_is_boundary = false;
        Walker hew = walker(v);
        if(boundary(*this, v))
        {
            if(val==2) return InvalidFaceID;
            vertex_is_boundary = true;
            hew = hew.circulate_vertex_ccw();
        }
        
        // Prepare some vectors for taking out the trash: We remove all old faces and all orphaned edges
        // and vertices
        vector<HalfEdgeID> garbage_halfedges;
        vector<FaceID> garbage_faces;
        vector<VertexID> garbage_vertices;
        
        vector<HalfEdgeID> loop; // The halfedges which form the outer loop of the merged one ring.
        
        // Below we loop over all faces adjacent to the vertex and add their halfedges to a big loop
        // which will form the loop of the new merged face. Below we remove from the loop edges
        // that appear twice (as each other's opposite).
        for(;!hew.full_circle(); hew = hew.circulate_vertex_ccw())
        {
            garbage_faces.push_back(hew.face());
            for(Walker hew2 = walker(hew.halfedge());
                !hew2.full_circle(); hew2 = hew2.circulate_face_ccw())
                loop.push_back(hew2.halfedge());
        }
        
        
        // Now merge the loops. We iteratively remove pairs of adjacent halfedges from the loop
        // if we find that the second is the opposite of the first since this is a degenerate
        // situation. However, we stop at four remaining halfedges since otherwise the loop degenerates
        // to zero after the next step if these four are also pairwise each other's opposites.
        int did_work;
        do
        {
            did_work = 0;
            vector<HalfEdgeID> loop_tmp(0);
            for(size_t i=0;i<loop.size();++i)
                if(walker(loop[i]).opp().halfedge() == loop[(i+1)%loop.size()])
                {
                    VertexID vid = walker(loop[i]).vertex();
                    if(vid != v)
                        garbage_vertices.push_back(walker(loop[i]).vertex());
                    garbage_halfedges.push_back(loop[i]);
                    garbage_halfedges.push_back(loop[(i+1)%loop.size()]);
                    loop[i] = InvalidHalfEdgeID;
                    loop[(i+1)%loop.size()] = InvalidHalfEdgeID;
                    ++did_work;
                    ++i;
                }
            for(size_t i=0;i<loop.size();++i)
                if(loop[i] != InvalidHalfEdgeID)
                    loop_tmp.push_back(loop[i]);
            loop = loop_tmp;
        } while(did_work > 0 && loop.size() > 4);
        
        // Check whether the loop is too long
        float loop_len=0.0;
        for(size_t i=0;i<loop.size();++i)
            loop_len += length(*this, loop[i]);
        if(loop_len > max_loop_length)
            return InvalidFaceID;
        
        // The following block checks wheteher the same halfedge appears twice. In this
        // case we fail since it means that the one ring contains the same face twice.
        vector<HalfEdgeID> loop_tmp = loop;
        sort(loop_tmp.begin(), loop_tmp.end());
        vector<HalfEdgeID>::iterator end_iter = unique(loop_tmp.begin(), loop_tmp.end());
        if(end_iter != loop_tmp.end())
            return InvalidFaceID;
        
        // Remove all faces and connected halfedges and the original vertex v.
        for(size_t i=0;i<garbage_vertices.size(); ++i)
            kernel.remove_vertex(garbage_vertices[i]);
        for(size_t i=0;i<garbage_faces.size(); ++i)
            kernel.remove_face(garbage_faces[i]);
        for(size_t i=0;i<garbage_halfedges.size(); ++i)
            kernel.remove_halfedge(garbage_halfedges[i]);
        if(!vertex_is_boundary)
            kernel.remove_vertex(v);
        
        // Create a new face for the merged one ring and link up all the halfedges
        // in the loop
        FaceID f = kernel.add_face();
        kernel.set_last(f,loop[0]);
        for(size_t i=0;i<loop.size(); ++i)
        {
            kernel.set_face(loop[i], f);
            Walker hw = walker(loop[i]);
            kernel.set_out(hw.vertex(),hw.opp().halfedge());
            link(loop[i],loop[(i+1)%loop.size()]);
            assert(hw.vertex() == walker(loop[(i+1)%loop.size()]).opp().vertex());
        }
        
        // Finally, we ensure boundary consitency for all vertices in the loop.
        for(size_t i=0;i<loop.size(); ++i)
        {
            Walker hw = walker(loop[i]);
            ensure_boundary_consistency(hw.vertex());
        }
        
        // Return the brand new merged face.
        return f;
    }
    
    
    
    bool Manifold::merge_faces(FaceID f, HalfEdgeID h)
    {
        //assert that we're merging a valid face with the corresponding halfedge
        assert(kernel.face(h) == f);
        
        HalfEdgeID ho = kernel.opp(h);
        FaceID fo = kernel.face(ho);
        HalfEdgeID hn = kernel.next(h);
        HalfEdgeID hon = kernel.next(ho);
        
        if(fo == f)
            return false;
        
        //boundary case
        if(kernel.vert(hn) == kernel.vert(hon))
            return false;
        
        HalfEdgeID hp = kernel.prev(h);
        HalfEdgeID hop = kernel.prev(ho);
        VertexID v = kernel.vert(h);
        VertexID vo = kernel.vert(ho);
        
        if(valency(*this, v) < 3 || valency(*this, vo) < 3)
            return false;
        
        link(hop, hn);
        link(hp, hon);
        kernel.set_out(vo, hon);
        kernel.set_out(v, hn);
        kernel.set_last(f, hn);
        
        HalfEdgeID hx = hon;
        
        assert(kernel.face(hx) == fo);
        while(kernel.face(hx) != f){
            kernel.set_face(hx, f);
            hx = kernel.next(hx);
        }
        
        ensure_boundary_consistency(v);
        ensure_boundary_consistency(vo);
        
        kernel.remove_halfedge(h);
        kernel.remove_halfedge(ho);
        kernel.remove_face(fo);
        
        return true;
    }
    
    FaceID Manifold::close_hole(HalfEdgeID h)
    {
        // invalid face is a hole
        if(kernel.face(h) == InvalidFaceID){
            FaceID f = kernel.add_face();
            kernel.set_last(f, h);
            do{
                kernel.set_face(h, f);
                h = kernel.next(h);
            }
            while(kernel.face(h) != f);
            return f;
        }
        return kernel.face(h);
    }
    
    VertexID Manifold::slit_vertex(VertexID v, HalfEdgeID h_in, HalfEdgeID h_out)
    {
        assert(kernel.face(h_in) != InvalidFaceID);
        assert(kernel.face(h_out) != InvalidFaceID);
        assert(kernel.opp(h_out) != h_in);
        
        VertexID v_new = kernel.add_vertex();
        pos(v_new) = pos(v);
        HalfEdgeID h = kernel.prev(h_out);
        kernel.set_vert(h, v_new);
        while ( h != h_in) {
            h = kernel.prev(kernel.opp(h));
            kernel.set_vert(h, v_new);
        }
        
        HalfEdgeID h_in_opp = kernel.opp(h_in);
        HalfEdgeID hn_in, hn_in_opp;
        if(kernel.face(h_in_opp) != InvalidFaceID)
        {
            hn_in = kernel.add_halfedge();
            kernel.set_face(hn_in, InvalidFaceID);
            glue(h_in_opp, hn_in);
            
            hn_in_opp = kernel.add_halfedge();
            kernel.set_face(hn_in_opp, InvalidFaceID);
            glue(h_in, hn_in_opp);
            
            link(hn_in_opp, hn_in);
            
            VertexID v_i = kernel.vert(h_in_opp);
            kernel.set_vert(hn_in_opp, v_i);
            kernel.set_out(v_i, hn_in);
        }
        else
        {
            hn_in_opp = h_in_opp;
            hn_in = kernel.prev(hn_in_opp);
            h_in_opp = kernel.opp(hn_in);
        }
        
        HalfEdgeID h_out_opp = kernel.opp(h_out);
        HalfEdgeID hn_out,hn_out_opp;
        if(kernel.face(h_out_opp) != InvalidFaceID)
        {
            hn_out = kernel.add_halfedge();
            kernel.set_face(hn_out, InvalidFaceID);
            glue(h_out_opp, hn_out);
            
            hn_out_opp = kernel.add_halfedge();
            kernel.set_face(hn_out_opp, InvalidFaceID);
            glue(h_out, hn_out_opp);

            link(hn_out, hn_out_opp);

            VertexID v_o = kernel.vert(h_out);
            kernel.set_vert(hn_out, v_o);
            kernel.set_out(v_o, hn_out_opp);
        }
        else
        {
            hn_out_opp = h_out_opp;
            hn_out = kernel.next(hn_out_opp);
            h_out_opp = kernel.opp(hn_out);
        }
        
        link(hn_out_opp, hn_in_opp);
        link(hn_in, hn_out);
        
        kernel.set_vert(hn_in, v);
        kernel.set_vert(hn_out_opp, v_new);
        
        kernel.set_out(v, hn_out);
        kernel.set_out(v_new, hn_in_opp);
        
        return v_new;
    }

    
    HalfEdgeID Manifold::slit_edges(VertexAttributeVector<int>& insel)
    {
        HalfEdgeID h;
        for(auto vid : vertices())
        {
            if(insel[vid])
            {
                HalfEdgeID h_in = InvalidHalfEdgeID, h_out = InvalidHalfEdgeID;
                Walker w = walker(vid);
                while(!w.full_circle())
                {
                    if(insel[w.vertex()]) {
                        if(h_in == InvalidHalfEdgeID) {
                            if(w.opp().face() == InvalidFaceID)
                                h_in = w.opp().next().opp().halfedge();
                            else
                                h_in = w.opp().halfedge();
                        }
                        else {
                            if(w.face() == InvalidFaceID)
                                h_out = w.prev().opp().halfedge();
                            else
                                h_out = w.halfedge();
                            break;
                        }
                    }
                    w = w.circulate_vertex_ccw();
                }
                if(h_in != InvalidHalfEdgeID &&
                   h_out != InvalidHalfEdgeID) {
                    VertexID v_new = slit_vertex(vid, h_in, h_out);
                    h = walker(v_new).halfedge();
                }
            }
        }
        return h;
    }

    
    void Manifold::flip_edge(HalfEdgeID h)
    {
        HalfEdgeID hn = kernel.next(h);
        HalfEdgeID hp = kernel.prev(h);
        HalfEdgeID ho = kernel.opp(h);
        HalfEdgeID hon = kernel.next(ho);
        HalfEdgeID hop = kernel.prev(ho);
        
        FaceID hf = kernel.face(h);
        FaceID hof = kernel.face(ho);
        
        VertexID hv = kernel.vert(h);
        VertexID hnv = kernel.vert(hn);
        VertexID hov = kernel.vert(ho);
        VertexID honv = kernel.vert(hon);
        
        // update face connectivity of halfedges changing face
        kernel.set_face(hop, hf);
        kernel.set_face(hp, hof);
        
        // connectivity of faces with halfedges of flipped edge
        kernel.set_last(hf, h);
        kernel.set_last(hof, ho);
        
        // connectivity of halfedges of first face after flip
        link(hn, h);
        link(h, hop);
        link(hop, hn);
        
        // connectivity of halfedges of second face after flip
        link(hon, ho);
        link(ho, hp);
        link(hp, hon);
        
        // connectivity of flipped edge and corresponding vertices
        kernel.set_vert(h, honv);
        kernel.set_vert(ho, hnv);
        
        if(kernel.out(hv) == ho)
            kernel.set_out(hv, hn);
        if(kernel.out(hov) == h)
            kernel.set_out(hov, hon);
        
        //
        //        // if the flip occurs next to a boundary, ensure the boundary consistency
        //        ensure_boundary_consistency(hv);
        //        ensure_boundary_consistency(hov);
    }
    
    
    /**********************************************
     * Private functions
     **********************************************/
    
    template<typename size_type, typename float_type, typename int_type>
    void Manifold::build_template(  size_type no_vertices,
                                  const float_type* vertvec,
                                  size_type no_faces,
                                  const int_type* facevec,
                                  const int_type* indices)
    {
        vector<VertexID> vids(no_vertices);
        
        // create vertices corresponding to positions stored in vertvec
        for(size_t i=0;i<no_vertices;++i)
        {
            const float_type* v0 = &vertvec[i*3];
            pos(vids[i] = kernel.add_vertex()) = Manifold::Vec(v0[0], v0[1], v0[2]);
        }
        
        auto hash_fun = [](const EdgeKey& k) {return k.hash();};
        //map over the created edges - needed to preserve edge uniqueness
		typedef unordered_map<EdgeKey, Edge, function<size_t(const EdgeKey&)>> EdgeMap;
        EdgeMap edge_map(no_vertices+no_faces,hash_fun);
        
        // counter that jumps between faces in indices vector
        int_type n  = 0;
        
        // create faces correspponding to faces stored in facevec
        for(size_type i = 0; i < no_faces; ++i){
            //amount of vertices in current face
            size_type N = facevec[i];
            vector<HalfEdgeID> fh;
            
            //each face indice results corresponds to 1 edge, 2 halfedges
            for(size_type j = 0; j < N; ++j){
                // ensure indice integrity
                
                assert(static_cast<size_type>(indices[j + n]) < no_vertices);
                assert(static_cast<size_type>(indices[(j + 1) % N + n]) < no_vertices);
                
                
                // each iteration uses two indices from the face
                VertexID v0 = vids[static_cast<size_type>(indices[j + n])];
                VertexID v1 = vids[static_cast<size_type>(indices[(j + 1) % N + n])];
                
                // create key and search map for edge
                EdgeKey ek(v0, v1);
                typename EdgeMap::iterator em_iter = edge_map.find(ek);
                
                // current edge has not been created
                if(em_iter == edge_map.end()){
                    // create edge for map
                    Edge e;
                    e.h0 = kernel.add_halfedge();
                    e.h1 = kernel.add_halfedge();
                    e.count = 1;
                    
                    // glue operation: 1 edge = 2 glued halfedges
                    glue(e.h0, e.h1);
                    
                    // update vertices with their outgoing halfedges
                    kernel.set_out(v0, e.h0);
                    kernel.set_out(v1, e.h1);
                    
                    // update halfedges with the vertices they point to
                    kernel.set_vert(e.h0, v1);
                    kernel.set_vert(e.h1, v0);
                    
                    // update map
                    edge_map[ek] = e;
                    
                    // update vector of halfedges belonging to current face
                    fh.push_back(e.h0);
                }
                else{
                    // update current face with halfedge from edge
                    fh.push_back(em_iter->second.h1);
                    
                    // asserting that a halfedge is visited exactly twice;
                    // once for each face on either side of the edge.
                    em_iter->second.count++;
                    assert(em_iter->second.count == 2);
                }
            }
            
            FaceID fid = kernel.add_face();
            for(size_type j = 0; j < N; ++j){
                // update halfedge with face
                kernel.set_face(fh[j], fid);
                
                // link operation: link two halfedges in the same face
                link(fh[j], fh[(j + 1) % N]);
            }
            //update face with the first halfedge created
            kernel.set_last(fid, fh[0]);
            
            // step to first index of next face
            n += N;
        }
        
        // test for unused vertices
        for(VertexIDIterator v = vertices_begin(); v != vertices_end(); ++v){
            assert( (*v) != InvalidVertexID);
            if(kernel.out(*v) == InvalidHalfEdgeID)
                kernel.remove_vertex(*v);
        }
        
        // boundary check while avoiding unused vertices
        for(VertexIDIterator v = vertices_begin(); v != vertices_end(); ++v){
            if((*v) != InvalidVertexID && in_use(*v))
                ensure_boundary_consistency(*v);
        }
    }
    void Manifold::link(HalfEdgeID h0, HalfEdgeID h1)
    {
        kernel.set_next(h0, h1);
        kernel.set_prev(h1, h0);
    }
    void Manifold::glue(HalfEdgeID h0, HalfEdgeID h1)
    {
        kernel.set_opp(h0, h1);
        kernel.set_opp(h1, h0);
    }
    void Manifold::ensure_boundary_consistency(VertexID v)
    {
        // boundary consistency check by performing two vertex circulations
        HalfEdgeID h = kernel.out(v);
        HalfEdgeID last = h;
        
        int c = 0;
        // step 1: circle through edges pointing away from vertex until reaching a null face
        while(kernel.face(h) != InvalidFaceID){
            h = kernel.opp(kernel.prev(h));
            if(h == last || ++c == 1e6) // We came full circle - vertex not boundary - return.
                return;
        }
        // null face encountered, we update our vertex with half edge index and prepare for step 2
        kernel.set_out(v, h);
        HalfEdgeID ho = kernel.opp(h);
        
        // step 2: circle through edges pointing towards vertex until reaching a null face
        while(kernel.face(ho) != InvalidFaceID){
            ho = kernel.opp(kernel.next(ho));
        }
        // null face encountered again, we update our edge with vertex index
        kernel.set_vert(ho, v);
        
        // remaining step is to make the in and out going edges link to each other
        link(ho, h);
    }
    void Manifold::remove_face_if_degenerate(HalfEdgeID h)
    {
        // face is degenerate if there is only two halfedges in face loop
        if(kernel.next(kernel.next(h)) == h)
        {
            HalfEdgeID hn = kernel.next(h);
            HalfEdgeID ho = kernel.opp(h);
            HalfEdgeID hno = kernel.opp(hn);
            VertexID hv = kernel.vert(h);
            VertexID hnv = kernel.vert(hn);
            FaceID f = kernel.face(h);
            
            assert(ho != hn);
            assert(h != hno);
            assert(hv != hnv);
            
            // glue opposite halfedge to halfedge opposite next halfedge, obsoleting h and hn from loop
            glue(ho, hno);
            
            // make v and vn point to opposite and next opposite halfedge, obsoleting h and hn from loop
            kernel.set_out(hnv, hno);
            kernel.set_out(hv, ho);
            
            // if face owning h is valid, remove face
            if(f != InvalidFaceID)
                kernel.remove_face(f);
            // remove the two invalid halfedges and the invalid face loop
            kernel.remove_halfedge(h);
            kernel.remove_halfedge(hn);
            
            // verify that v and vn is sane after removing the degenerate face
            ensure_boundary_consistency(hv);
            ensure_boundary_consistency(hnv);
        }
    }
    
    vector<HalfEdgeID> Manifold::bridge_faces(FaceID f0, FaceID f1, const vector<pair<VertexID, VertexID> >& pairs)
    {
        // Let N be the number of vertex pairs to connect by edges
        size_t N = pairs.size();
        
        // We now create N pairs of edges, glue each pair and let them point to the appropriate
        // vertices.
        vector<HalfEdgeID> new_halfedges(N);
        vector<HalfEdgeID> new_halfedges_opp(N);
        for(size_t i=0;i<N; ++i)
        {
            new_halfedges[i] = kernel.add_halfedge();
            new_halfedges_opp[i] = kernel.add_halfedge();
            glue(new_halfedges[i], new_halfedges_opp[i]);
            kernel.set_vert(new_halfedges[i], pairs[i].second);
            kernel.set_vert(new_halfedges_opp[i], pairs[i].first);
        }
        
        // We now maintain some halfedge indices to right before
        // and right after the point we are trying to connect on
        // each of the two loops.
        HalfEdgeID h0p = kernel.last(f0);
        HalfEdgeID h1n = kernel.last(f1);
        
        // loop over all connections and link the new halfedges with the old
        // ones.
        for(size_t i=0;i<N; ++i)
        {
            while(kernel.vert(h0p) != pairs[i].first)
                h0p = kernel.next(h0p);
            while(kernel.vert(kernel.prev(h1n)) != pairs[i].second)
                h1n = kernel.prev(h1n);
            
            HalfEdgeID h0n = kernel.next(h0p);
            HalfEdgeID h1p = kernel.prev(h1n);
            
            link(h0p, new_halfedges[i]);
            link(new_halfedges[i],h1n);
            link(h1p, new_halfedges_opp[i]);
            link(new_halfedges_opp[i],h0n);
            
            h0p = new_halfedges_opp[i];
            h1n = new_halfedges_opp[i];
        }
        
        // Create the faces and their connections
        for(size_t i=0;i<N; ++i)
        {
            HalfEdgeID last = new_halfedges[i];
            FaceID f = kernel.add_face();
            kernel.set_last(f, last);
            HalfEdgeID h = last;
            do
            {
                kernel.set_face(h, f);
                h = kernel.next(h);
            } while(h != last);
        }
        
        // Delete the old faces that were bridged.
        kernel.remove_face(f0);
        kernel.remove_face(f1);
        
        new_halfedges.insert(new_halfedges.end(), new_halfedges_opp.begin(), new_halfedges_opp.end());
        return new_halfedges;
    }
    
    /***************************************************
     * Namespace functions
     ***************************************************/
    bool valid(const Manifold& m)
    {
        // Verify components of halfedges
        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h){
            Walker j = m.walker(*h);
            
            if(j.vertex() == InvalidVertexID){
                cout << "Halfedge lacks vert" << endl;
                return false;
            }
            if(j.next().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks next" << endl;
                return false;
            }
            if(j.prev().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks prev" << endl;
                return false;
            }
            if(j.opp().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks opp" << endl;
                return false;
            }
            
        }
        // Verify components of vertices
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            vector<VertexID> link;
            
            // circulate the halfedges of vertex
            for(Walker j = m.walker(*v); !j.full_circle(); j = j.circulate_vertex_cw()){
                // test halfedges around v
                if(j.halfedge() == InvalidHalfEdgeID){
                    cout << "Vertex circulation produced invalid halfedge" << endl;
                    return false;
                }
                VertexID ring_v = j.vertex();
                
                // test one-ring for multiple occurences of vertex
                if(find(link.begin(), link.end(), ring_v) != link.end()){
                    cout << "Vertex appears two times in one-ring of vertex" << endl;
                    return false;
                }
                link.push_back(ring_v);
                
                // test for infinite loop around vertex
                if(static_cast<size_t>(j.no_steps()) > m.no_vertices()){
                    cout << "Vertex loop CW contains more vertices than manifold" << endl;
                    return false;
                }
            }
            
            for(Walker j = m.walker(*v); !j.full_circle(); j = j.circulate_vertex_ccw()) {
                if(static_cast<size_t>(j.no_steps()) > m.no_vertices()){
                    cout << "Vertex loop CCW contains more vertices than manifold" << endl;
                    return false;
                }
            }
            
            // test one_ring size for boundary consistency
            if(link.size() <= 2){
                Walker j = m.walker(*v);
                
                if(j.face() != InvalidFaceID && j.opp().face() != InvalidFaceID)
                {
                    if(link.size()==1)
                        cout << "Vertex contains only a single incident edge" << endl;
                    //
                    //                    cout << "Vertex contains only " << link.size() <<" incident edges" << endl;
                }
                else
                    cout << "Boundary vertex contains only " << link.size() <<" incident edges" << endl;
            }
        }
        // verify components of faces
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            // count edges on face
            Walker j = m.walker(*f);
            
            for(; !j.full_circle(); j = j.circulate_face_cw()){
                // test that all halfedges in faces bind properly to their face
                if(j.face() != *f){
                    cout << "Face is inconsistent, halfedge is not bound to face" << endl;
                    return false;
                }
            }
            // test faces for valid geometrical properties
            if(j.no_steps() < 3){
                cout << "Face contains less than 3 edges" << endl;
                return false;
            }
            // test for infinite loop around face
            if(j.no_steps() > m.no_halfedges() * 0.5f){
                cout << "Face loop contains more halfedges than manifold" << endl;
                return false;
            }
        }
        return true;
    }
    
    void bbox(const Manifold& m, Manifold::Vec& pmin, Manifold::Vec& pmax)
    {
        if(m.no_vertices()==0)
            return;
        VertexIDIterator v = m.vertices_begin();
        pmin = pmax = m.pos(*v);
        ++v;
        for(; v != m.vertices_end(); ++v){
            pmin = v_min(m.pos(*v), pmin);
            pmax = v_max(m.pos(*v), pmax);
        }
    }
    
    void bsphere(const Manifold& m, Manifold::Vec& c, float& r)
    {
        Manifold::Vec p0, p7;
        bbox(m, p0, p7);
        Manifold::Vec rad = (p7 - p0) * 0.5f;
        c = p0 + rad;
        r = rad.length();
    }
    
    
    
    bool precond_collapse_edge(const Manifold& m, HalfEdgeID h)
    {
        Walker hew = m.walker(h);
        HalfEdgeID ho = hew.opp().halfedge();
        VertexID v0 = hew.opp().vertex();
        VertexID v1 = hew.vertex();
        
        // get the one-ring of v0
        vector<VertexID> link0;
        vector<FaceID> faces0;
        int k = 0;
        for(Walker vj = m.walker(h);
            !vj.full_circle(); vj = vj.circulate_vertex_ccw()){
            link0.push_back(vj.vertex());
            if(vj.halfedge() != h)
                faces0.push_back(vj.face());
            if(++k>m.no_vertices())
            {
                cout << "mesh is corrupted" << endl;
                return false;
            }
        }
        assert(link0.size() > 1);
        
        // get the one-ring of v1
        vector<VertexID> link1;
        vector<FaceID> faces1;
        k=0;
        for(Walker vj = m.walker(ho);
            !vj.full_circle(); vj = vj.circulate_vertex_ccw()){
            link1.push_back(vj.vertex());
            if(vj.halfedge() != ho)
                faces1.push_back(vj.face());
            if(++k>m.no_vertices())
            {
                cout << "mesh is corrupted" << endl;
                return false;
            }
        }
        assert(link1.size() > 1);
        
        // sort the vertices of the two rings
        sort(link0.begin(), link0.end());
        sort(link1.begin(), link1.end());
        
        // get the intersection of the shared vertices from both rings
        vector<VertexID> lisect;
        std::back_insert_iterator<vector<VertexID> > lii(lisect);
        set_intersection(link0.begin(), link0.end(),
                         link1.begin(), link1.end(),
                         lii);
        
        // sort the faces of the two rings
        sort(faces0.begin(), faces0.end());
        sort(faces1.begin(), faces1.end());
        
        // get the intersection of the shared faces from both rings
        vector<FaceID> fisect;
        std::back_insert_iterator<vector<FaceID> > fii(fisect);
        set_intersection(faces0.begin(), faces0.end(),
                         faces1.begin(), faces1.end(),
                         fii);
        if(fisect.size() > 0)
            return false;
        
         k = 0;
        // if the adjacent face is a triangle (see 2)
        if(hew.next().next().next().halfedge() == h){
            VertexID v = hew.next().vertex();
            
            // valency test (see 5)
            if(valency(m, v) < 3)
                return false;
            
            // remove the vertex shared by the two rings from the intersection
            vector<VertexID>::iterator iter;
            iter = find(lisect.begin(), lisect.end(), v);
            assert(iter != lisect.end());
            lisect.erase(iter);
            ++k;
        }
        // if the adjacent face is a triangle (see 2)
        if(hew.opp().next().next().next().halfedge() == hew.opp().halfedge()){
            VertexID v = hew.opp().next().vertex();
            
            // valency test (see 5)
            if(valency(m, v) < 3)
                return false;
            
            // remove the vertex shared by the two rings from the intersection
            vector<VertexID>::iterator iter;
            iter = find(lisect.begin(), lisect.end(), v);
            assert(iter != lisect.end());
            lisect.erase(iter);
            ++k;
        }
        // double edge test (see 3)
        if(lisect.size() != 0)
            return false;
        
        // tetrahedon test (see 4)
        if(k == 2 && (link0.size() + link1.size() == 6))
            return false;
        
        // test that we do not merge holes (see 6)
        if(boundary(m, v0) && boundary(m, v1) && hew.face() != InvalidFaceID && hew.opp().face() != InvalidFaceID)
            return false;
        
        return true;
    }
    
    bool precond_flip_edge(const Manifold& m, HalfEdgeID h)
    {
        Walker j = m.walker(h);
        
        FaceID hf = j.face();
        FaceID hof = j.opp().face();
        
        // boundary case
        if(hf == InvalidFaceID || hof == InvalidFaceID)
           return false;
        
        
        // We can only flip an edge if both incident polygons are triangles.
        if(no_edges(m, hf) != 3 || no_edges(m, hof) !=3)
            return false;
        

        // non boundary vertices with a valency of less than 4(less than 3 after operation) degenerates mesh.
        VertexID hv = j.vertex();
        VertexID hov = j.opp().vertex();
        if((valency(m, hv) < 4 && !boundary(m, hv)) || (valency(m, hov) < 4 && !boundary(m, hov))){
            return false;
        }
        
        // Disallow flip if vertices being connected already are.
        VertexID hnv = j.next().vertex();
        VertexID honv = j.opp().next().vertex();
        if(connected(m, hnv, honv)){
           return false;
        }
        
        return true;
    }
    
    bool boundary(const Manifold& m, VertexID v)
    {
        Walker j  = m.walker(v);
        return boundary(m, j.halfedge());
    }

    int valency(const Manifold& m, VertexID v)
    {
        return circulate_vertex_ccw(m,v, (std::function<void(Walker&)>)[](Walker){});
    }
    
    Manifold::Vec normal(const Manifold& m, VertexID v)
    {
        Manifold::Vec p0 = m.pos(v);
        vector<Manifold::Vec> one_ring;
        
        // run through outgoing edges, and store them normalized
        circulate_vertex_ccw(m, v, (std::function<void(VertexID)>)[&](VertexID vn) {
            Manifold::Vec edge = m.pos(vn) - p0;
            double l = length(edge);
            if(l > 0.0)
                one_ring.push_back(edge/l);
        });
        int N = one_ring.size();
        if(N<2)
            return Manifold::Vec(0);
        
        size_t N_count = N;
        size_t N_start = 0;
        if(boundary(m, v))
            N_start = 1;
        
        // sum up the normals of each face surrounding the vertex
        Manifold::Vec n(0);
        for(size_t i = N_start; i < N_count; ++i){
            Manifold::Vec e0 = one_ring[i];
            Manifold::Vec e1 = one_ring[(i+1) % N];
            
            Manifold::Vec n_part = normalize(cross(e0, e1));
            n += n_part * acos(max(-1.0, min(1.0, dot(e0, e1))));
        }
        
        // normalize and return the normal
        float sqr_l = sqr_length(n);
        if(sqr_l > 0.0f) return n / sqrt(sqr_l);
        
        return n;
    }
    
    
    bool connected(const Manifold& m, VertexID v0, VertexID v1)
    {
        bool c=false;
        circulate_vertex_ccw(m, v0, (std::function<void(VertexID)>)[&](VertexID v){ c |= (v==v1);});
        return c;
    }

    
    int no_edges(const Manifold& m, FaceID f)
    {
        return circulate_face_ccw(m, f, (std::function<void(Walker&)>)[](Walker w){});
    }
    
    Manifold::Vec normal(const Manifold& m, FaceID f)
    {
        vector<Manifold::Vec> v;
        
        int k= circulate_face_ccw(m, f, (std::function<void(VertexID)>)[&](VertexID vid) {
            v.push_back(m.pos(vid));
        });
        
        Manifold::Vec norm(0);
        for(int i=0;i<k;++i)
        {
            norm[0] += (v[i][1]-v[(i+1)%k][1])*(v[i][2]+v[(i+1)%k][2]);
            norm[1] += (v[i][2]-v[(i+1)%k][2])*(v[i][0]+v[(i+1)%k][0]);
            norm[2] += (v[i][0]-v[(i+1)%k][0])*(v[i][1]+v[(i+1)%k][1]);
        }
        float l = sqr_length(norm);
        if(l>0.0f)
            norm /= sqrt(l);
        return norm;
    }
    
    
    double area(const Manifold& m, FaceID fid)
    {
        // Get all projected vertices
        vector<Manifold::Vec> vertices;
        int N = circulate_face_ccw(m, fid, (std::function<void(VertexID)>)[&](VertexID vid) {
            vertices.push_back(m.pos(vid));
        });

        
        double area = 0;
        Manifold::Vec norm = normal(m,fid);
        for(int i = 1; i < N-1; ++i)
            area += 0.5 * dot(norm,cross(vertices[i]-vertices[0], vertices[(i+1 )]-vertices[0]));
        return area;
    }
    
    Manifold::Vec centre(const Manifold& m, FaceID f)
    {
        Manifold::Vec c(0);
        int n = circulate_face_ccw(m, f, (std::function<void(VertexID)>)[&](VertexID v) {c+=m.pos(v);});
        return c / n;
    }
    
    double perimeter(const Manifold& m, FaceID f)
    {
        double l=0.0;
        circulate_face_ccw(m, f, (std::function<void(HalfEdgeID)>)[&](HalfEdgeID h) { l+= length(m, h);});
        return l;
    }
    
    bool boundary(const Manifold& m, HalfEdgeID h)
    {
        Walker w = m.walker(h);
        return w.face() == InvalidFaceID || w.opp().face() == InvalidFaceID;
    }
    
    double length(const Manifold& m, HalfEdgeID h)
    {
        Walker w = m.walker(h);
        return (m.pos(w.vertex()) - m.pos(w.opp().vertex())).length();
    }
}
