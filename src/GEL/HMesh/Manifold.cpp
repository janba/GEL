 /* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */
#include <GEL/Geometry/TriMesh.h>
#include <GEL/Geometry/bounding_sphere.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/cleanup.h>
#include <GEL/HMesh/triangulate.h>
#include <GEL/Util/Assert.h>

#include <GEL/Util/AssociativeContainers.h>

#include <iostream>
#include <iterator>
#include <map>
#include <ranges>
#include <stack>
#include <vector>

namespace HMesh
{
    using namespace std;
    using namespace Geometry;
    using namespace CGLA;

    /*********************************************
	 * Public functions
	 *********************************************/

    bool Manifold::remove_face(FaceID fid)
    {
        if(!in_use(fid))
            return false;
        
        HalfEdgeID h = kernel.last(fid);
        HalfEdgeID last = h;
        vector<HalfEdgeID> halfedges;
        do {
            halfedges.push_back(h);
            kernel.set_face(h, InvalidFaceID);
            h = kernel.next(h);
        }
        while(h != last);
        
        vector<HalfEdgeID> halfedges_garbage;
        vector<VertexID> vertices_garbage;
        for(size_t i=0;i<halfedges.size(); ++i)
        {
            HalfEdgeID h = halfedges[i];
            HalfEdgeID ho = kernel.opp(h);
            HalfEdgeID hn = kernel.next(h);
            HalfEdgeID hno = kernel.opp(hn);
            bool f_ho_invalid = kernel.face(ho) == InvalidFaceID;
            bool f_hno_invalid = kernel.face(hno) == InvalidFaceID;

            if(f_ho_invalid) {
                halfedges_garbage.push_back(h);
                halfedges_garbage.push_back(ho);
            }
            
            VertexID v1 = kernel.vert(h);
            if (f_ho_invalid || f_hno_invalid) {
                if (valency(v1) > 2) {
                    HalfEdgeID h_v1in = f_ho_invalid ? kernel.prev(ho) : h;
                    HalfEdgeID h_v1out = f_hno_invalid ? kernel.next(hno) : hn;
                    link(h_v1in, h_v1out);
                    kernel.set_out(v1, h_v1out);
                }
                else vertices_garbage.push_back(v1);
            }
        }

        for(size_t i=0;i<halfedges_garbage.size();++i){
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
        circulate_vertex_ccw(*this, vid, [&](FaceID f) {
            if(in_use(f))
                faces.push_back(f);
        });
        for(auto f: faces)
            remove_face(f);

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
			GEL_ASSERT_EQ(kernel.vert(kernel.opp(he)), hov);
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
    }

    FaceID Manifold::split_face_by_edge(FaceID f, VertexID v0, VertexID v1)
    {
        // make sure we are not trying to connect a vertex to itself
        GEL_ASSERT_NEQ(v0, v1);
        GEL_ASSERT_NEQ(f, InvalidFaceID);
        GEL_ASSERT_FALSE(connected(v0, v1));

        // Find the edge along f that emanates from v0
        HalfEdgeID h0 = kernel.out(v0);
        Walker w = walker(v0);
        for (int steps_v = 0; !w.full_circle(); w = w.circulate_vertex_ccw(), ++steps_v) {
            if (w.face() == f) {
                h0 = w.halfedge();
                break;
            }
        }
        // If v0 is not a corner of f, we bail out
        if (kernel.face(h0) != f) { return InvalidFaceID; }

        // the halfedge belonging to f, going out from v0, is denoted h. Move along h until we hit v1.
        // now we have the halfedge which belongs to f and points to v1.
        HalfEdgeID h = h0;
        while (kernel.vert(h) != v1) {
            h = kernel.next(h);
            if (h == h0) { // Oops came full circle, bail out...
                return InvalidFaceID;
            }
        }

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
        while (kernel.vert(h) != v0) {
            kernel.set_face(h, f2);
            h = kernel.next(h);
        }
        kernel.set_face(h, f2);
        GEL_ASSERT_NEQ(h, h1);

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
        GEL_ASSERT_EQ(kernel.next(kernel.opp(kernel.prev(h1))), h0);
        GEL_ASSERT_EQ(kernel.next(kernel.opp(kernel.prev(h0))), h1);
        GEL_ASSERT_EQ(kernel.face(kernel.next(hb)), f2);
        GEL_ASSERT_EQ(kernel.face(kernel.next(kernel.next(hb))), f2);
        GEL_ASSERT_EQ(kernel.face(hb), f2);

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

        // create the new vertex with middle of edge as position and update connectivity
        VertexID vn = kernel.add_vertex();
        positions[vn] = .5f * (positions[v] + positions[vo]);
        kernel.set_out(vn, h);

        // create two necessary halfedges, and update connectivity
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
        if (kernel.face(h) != InvalidFaceID) kernel.set_last(kernel.face(h), hn);

        if (kernel.face(ho) != InvalidFaceID) kernel.set_last(kernel.face(ho), ho);

        // update new edge with faces from existing edge
        kernel.set_face(hn, kernel.face(h));
        kernel.set_face(hno, kernel.face(ho));

        return vn;
    }

    FaceID Manifold::add_boundary_face(const HalfEdgeID h)
    {
        GEL_ASSERT(kernel.in_use(h), "h must be valid");
        GEL_ASSERT(kernel.face(h) == InvalidFaceID, "h must be on boundary");

        const FaceID f = add_face({Vec3d(0), Vec3d(0), Vec3d(0)});
        const auto hn0_o = kernel.opp(kernel.last(f));
        const auto hn1_o = kernel.next(hn0_o);
        // const auto hn2_o = kernel.next(hn1_o);
        GEL_ASSERT(stitch_boundary_edges(h, hn1_o));

        // face | kernel.last | kernel.vert points to the dangling vertex
        return f;
    }

    FaceID Manifold::add_boundary_face(const HalfEdgeID h0, const HalfEdgeID h1)
    {
        GEL_ASSERT(kernel.in_use(h0), "h0 must be valid");
        GEL_ASSERT(kernel.in_use(h1), "h0 must be valid");
        GEL_ASSERT(kernel.face(h0) == InvalidFaceID, "h0 must be on boundary");
        GEL_ASSERT(kernel.face(h1) == InvalidFaceID, "h1 must be on boundary");
        GEL_ASSERT(kernel.next(h0) == h1, "h1 must follow h0");

        const FaceID f = add_face({Vec3d(0), Vec3d(0), Vec3d(0)});
        const auto hn0_o = kernel.opp(kernel.last(f));
        const auto hn1_o = kernel.next(hn0_o);
        const auto hn2_o = kernel.next(hn1_o);

        GEL_ASSERT(stitch_boundary_edges(h0, hn2_o));
        GEL_ASSERT(stitch_boundary_edges(h1, hn1_o));

        // face | kernel.last | kernel.opp is the new boundary halfedge
        return f;
    }

    VertexID Manifold::split_vertex(HalfEdgeID h_in, HalfEdgeID h_out)
    {
        GEL_ASSERT(kernel.in_use(h_in), "h1 must be valid");
        GEL_ASSERT(kernel.in_use(h_in), "h2 must be valid");
        GEL_ASSERT_EQ(kernel.vert(kernel.opp(h_out)), kernel.vert(h_in),
                      "h_in must point to h_out's origin");
        const auto v = kernel.vert(h_in);

        // If the halfedges are on the border, then we don't have to perform a slit, we can just directly add a face.
        if (kernel.face(h_in) == InvalidFaceID && kernel.face(h_out) == InvalidFaceID) {
            const auto f = add_boundary_face(h_in);
            const auto vn = kernel.vert(kernel.last(f));
            pos(vn) = pos(v);
            const auto hn = kernel.opp(kernel.last(f));
            add_boundary_face(hn, h_out);

            return vn;
        } else if (kernel.face(h_out) == InvalidFaceID) {
            auto h_out_opp = kernel.opp(h_out);
            auto f_dummy = add_boundary_face(kernel.prev(h_out), h_out);

            GEL_ASSERT_NEQ(kernel.face(h_in), InvalidFaceID);
            GEL_ASSERT_NEQ(kernel.face(h_out_opp), InvalidFaceID);

            const auto vn = slit_vertex_impl(h_in, kernel.opp(h_out_opp));
            GEL_ASSERT_NEQ(vn, InvalidVertexID);

            const auto f = close_hole(kernel.out(vn));

            split_face_by_edge(f, v, vn);

            GEL_ASSERT(remove_face(f_dummy));

            return vn;
        } else if (kernel.face(h_in) == InvalidFaceID) {
            auto h_in_opp = kernel.opp(h_in);
            auto f_dummy = add_boundary_face(h_in, kernel.next(h_in));

            GEL_ASSERT_NEQ(kernel.face(h_in_opp), InvalidFaceID);
            GEL_ASSERT_NEQ(kernel.face(h_out), InvalidFaceID);

            const auto vn = slit_vertex_impl(kernel.opp(h_in_opp), h_out);
            GEL_ASSERT_NEQ(vn, InvalidVertexID);

            const auto f = close_hole(kernel.out(vn));

            split_face_by_edge(f, v, vn);

            GEL_ASSERT(remove_face(f_dummy));

            return vn;
        } else if ((kernel.face(kernel.opp(h_in)) == InvalidFaceID) != (kernel.face(kernel.opp(h_out)) == InvalidFaceID)) {
            const auto v_new = slit_vertex_impl(h_in, h_out);
            GEL_ASSERT_NEQ(v_new, InvalidVertexID);

            if (kernel.vert(kernel.next(kernel.opp(h_in))) == v) {
                const auto h_in_opp = kernel.opp(h_in);
                const auto h_in_opp_next = kernel.next(h_in_opp);
                const auto h_in_opp_prev = kernel.prev(h_in_opp);
                const auto f = add_boundary_face(h_in_opp, h_in_opp_next);
                const auto h_b = kernel.opp(kernel.last(f));
                add_boundary_face(h_in_opp_prev, h_b);
            } else {
                const auto h_out_opp = kernel.opp(h_out);
                const auto h_out_opp_prev = kernel.prev(h_out_opp);
                const auto h_out_opp_next = kernel.next(h_out_opp);
                const auto f = add_boundary_face(h_out_opp_prev, h_out_opp);
                const auto h_b = kernel.opp(kernel.last(f));
                add_boundary_face(h_b, h_out_opp_next);
            }
            return v_new;
        } else {
            const auto vn = slit_vertex_impl(h_in, h_out);
            GEL_ASSERT_NEQ(vn, InvalidVertexID);

            GEL_ASSERT_NEQ(kernel.face(h_in), InvalidFaceID);
            GEL_ASSERT_NEQ(kernel.face(h_out), InvalidFaceID);

            const auto f = close_hole(kernel.out(vn));

            split_face_by_edge(f, v, vn);

            return vn;
        }
    }

    VertexSet link_intersection(const Manifold& m, VertexID v0, VertexID v1)
    {
        // get the one-ring of v0
        vector<VertexID> link0;
        circulate_vertex_ccw(m, v0, [&](VertexID vn) {
            link0.emplace_back(vn);
        });
		
        // get the one-ring of v1
        vector<VertexID> link1;
        circulate_vertex_ccw(m, v1, [&](VertexID vn) {
            link1.emplace_back(vn);
        });
		
        // sort the vertices of the two rings
        sort(link0.begin(), link0.end());
        sort(link1.begin(), link1.end());
		
        // get the intersection of the shared vertices from both rings
        
        VertexSet vset;
        set_intersection(link0.begin(), link0.end(),
						 link1.begin(), link1.end(),
						 inserter(vset,vset.begin()));
        
        return vset;
    }
    
    bool Manifold::stitch_boundary_edges(HalfEdgeID h0, HalfEdgeID h1)
    {
        GEL_ASSERT(kernel.in_use(h0), "h0 must be valid");
        GEL_ASSERT(kernel.in_use(h1), "h1 must be valid");
        // Cannot stitch an edge with itself
        if(h0 == h1) {
        //    cout << "Stupid, you are stitching edge with self!" << endl;
            return false;
        }
        // Only stitch a pair of boundary edges.
        if(kernel.face(h0) == InvalidFaceID && kernel.face(h1) == InvalidFaceID)
        {
            HalfEdgeID h0o = kernel.opp(h0);
            HalfEdgeID h1o = kernel.opp(h1);
            VertexID v0a = kernel.vert(h0);
            VertexID v0b = kernel.vert(kernel.opp(h1));
            VertexID v1b = kernel.vert(h1);
            VertexID v1a = kernel.vert(kernel.opp(h0));
            FaceID f0 = kernel.face(h0o);
            FaceID f1 = kernel.face(h1o);

            // We refuse to stitch two boundary edges belonging to the same face.
            if (f0 == f1) {
            //    cout << "Douche it is the same face!" << endl;
                return false;
            }

            //If the vertices are already connected, welding them together will be awkward.
            if(connected(v0a, v0b) || connected(v1a, v1b)) {
            //    cout << "Stupid, end points are connected!" << endl;
                return false;
            }
            if((v0a != v0b && v1a != v1b) &&
               (connected(v0a, v1b) || connected(v1a, v0b))) {
            //    cout << "Stupid, end points are distinct and cross connected!" << endl;
                return false;
            }
            
            // This check ensures that we do not create a valency 1 vertex by welding
            // two edges whose opposite edges have the property that second one is the next
            // edge of the first one or vice versa.
            if(v1a == v1b && kernel.next(h0o) == h1o) {
            //    cout << "Silly: would have made val 1" << endl;
                return false;
            }
            if(v0a == v0b && kernel.next(h1o) == h0o) {
            //    cout << "WE CANT HAVE VALENCE 1!" << endl;
                return false;
            }

            // If we are here, we can stitch the edges, but it may be the case that
            // the two edges share either or both vertices.
            // For each endpoint X=0 or X=1,  we check whether vXa and vXb are distinct.
            // If they are distinct, we ensure that the incoming halfedges of vXb point to vXa.
            // We also link the two boundary edges h1p and h0n since h1p is now the successor to
            // h0n. Likewise, for h0p and h1n. Finally, we remove the spurious vertex vXb.
            if(v0a != v0b)
            {
                circulate_vertex_ccw(*this, v0b, [&](Walker& hew) {
                    kernel.set_vert(hew.opp().halfedge(), v0a);
                });

                HalfEdgeID h1p = kernel.prev(h1);
                HalfEdgeID h0n = kernel.next(h0);
                link(h1p, h0n);
                kernel.remove_vertex(v0b);
            }
            
            if(v1a != v1b)
            {
                circulate_vertex_ccw(*this, v1b, [&](Walker& hew) {
                    kernel.set_vert(hew.opp().halfedge(), v1a);
                });

                HalfEdgeID h0p = kernel.prev(h0);
                HalfEdgeID h1n = kernel.next(h1);
                link(h0p, h1n);
                kernel.remove_vertex(v1b);
            }

            // Finally, we perform the actual gluing of the edges. Note that since there are always
            // edges along a boundary, it is actually the two opposite edges that are being glued, and the
            // two input edges are removed. The outgoing halfedges are updated for the two vertices at
            // either end of the glued halfedges.
            glue(h0o, h1o);
            kernel.remove_halfedge(h0);
            kernel.remove_halfedge(h1);
            kernel.set_out(v1a, h1o);
            kernel.set_out(v0a, h0o);
            
            return true;
        }
    //    cout << "Idiot! Those are not both boundary edges" << endl;
        return false;
    }
    
    bool Manifold::merge_boundary_vertices(VertexID v0, VertexID v1) {
        auto find_hi_ho = [this](VertexID v, HalfEdgeID& hi, HalfEdgeID& ho) -> bool {
            int sanity_count = 0;
            circulate_vertex_ccw(*this, v, [&](Walker w) {
                if(w.face() == InvalidFaceID) {
                    ho = w.halfedge();
                    sanity_count += 1;
                }
                if(w.opp().face() == InvalidFaceID)
                    hi = w.opp().halfedge();
            });
            return sanity_count == 1;
        };
        HalfEdgeID h0i, h0o, h1i, h1o;
        
        vector<VertexID> r0;
        circulate_vertex_ccw(*this, v0, [&](VertexID v){r0.push_back(v); cout  << v << ",";});
        vector<VertexID> r1;
        circulate_vertex_ccw(*this, v1, [&](VertexID v){r1.push_back(v); cout  << v << ",";});

        if(find(begin(r0),end(r0),v1) != end(r0)) {
        //    cout << "Oops " << v1  << " in 1-ring of " << v0;
            return false;
        }
        if(find(begin(r1),end(r1),v0) != end(r1)) {
        //    cout << "Oops " << v0  << " in 1-ring of " << v1;
            return false;
        }
        
        sort(begin(r0), end(r0));
        sort(begin(r1), end(r1));
        vector<VertexID> risect;
        set_intersection(begin(r0), end(r0), begin(r1), end(r1), back_inserter(risect));
        if(!risect.empty())
        {
        //    cout << "One rings overlap" << endl;
            return false;
        }

        if (find_hi_ho(v0, h0i, h0o) && find_hi_ho(v1, h1i, h1o)){
            link(h0i, h1o);
            link(h1i, h0o);
            kernel.set_vert(h1i, v0);
            kernel.set_vert(kernel.opp(h1o), v0);
            HalfEdgeID h0 = kernel.opp(kernel.out(v0));
            HalfEdgeID h = h0;
            do {
                kernel.set_vert(h, v0);
                h = kernel.opp(kernel.next(h));
            }
            while( h != h0);
            kernel.remove_vertex(v1);
//            cout << "MERGING " << v0 << " and " << v1 << ", halfedge: " << h0i << "," << h0o << "," << h1i << "," << h1o << endl;
            return true;
        }
        cout << "Cannot merge " << v0 << " and " << v1 << endl;
        return false;
    }
    
    FaceID Manifold::merge_one_ring(VertexID v)
    {
        // If the vertex is either not in use or has just
        // one incident edge (or less), bail out.
        int val = valency(v);
        if(!in_use(v) || val<2)
            return InvalidFaceID;
        
        // If the vertex is  a boundary vertex, we preparte the walker so that the
        // first face visited is not the invalid face outside the boundary. If the boundary
        // vertex is adjacent to only one vertex, there is little to do and we bail out.
        bool vertex_is_boundary = false;
        HalfEdgeID h = boundary_edge(v);
        Walker hew = walker(v);
        if(h != InvalidHalfEdgeID)
        {
            if(val==2)
                return InvalidFaceID;
            vertex_is_boundary = true;
            hew = walker(h).circulate_vertex_ccw();
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
            Walker w = walker(hew.halfedge());
            for(; !w.full_circle(); w = w.next())
                loop.push_back(w.halfedge());
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
                    VertexID vid = kernel.vert(loop[i]);
                    if(vid != v)
                        garbage_vertices.push_back(vid);
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
        
        
        // The following block checks wheteher the same halfedge appears twice. In this
        // case we fail since it means that the one ring contains the same face twice.
        vector<HalfEdgeID> loop_tmp = loop;
        sort(loop_tmp.begin(), loop_tmp.end());
        vector<HalfEdgeID>::iterator end_iter = unique(loop_tmp.begin(), loop_tmp.end());
        if(end_iter != loop_tmp.end())
            return InvalidFaceID;
        
        // Remove all faces and connected halfedges and the original vertex v.
        for(auto v: garbage_vertices)
            kernel.remove_vertex(v);
        for(auto f: garbage_faces)
            kernel.remove_face(f);
        for(auto h: garbage_halfedges)
            kernel.remove_halfedge(h);
        if(!vertex_is_boundary)
            kernel.remove_vertex(v);
        
        // Create a new face for the merged one ring and link up all the halfedges
        // in the loop
        FaceID f = kernel.add_face();
        kernel.set_last(f,loop[0]);
        for(size_t i=0;i<loop.size(); ++i)
        {
            kernel.set_face(loop[i], f);
            kernel.set_out(kernel.vert(loop[i]),kernel.opp(loop[i]));
            link(loop[i],loop[(i+1)%loop.size()]);
        }
        
        // Return the brand new merged face.
        return f;
    }
    
    
    
    bool Manifold::merge_faces(FaceID f, HalfEdgeID h)
    {
        //assert that we're merging a valid face with the corresponding halfedge
        GEL_ASSERT_EQ(kernel.face(h), f);
        
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
        
        if(valency(v) < 3 || valency(vo) < 3)
            return false;
        
        link(hop, hn);
        link(hp, hon);
        kernel.set_out(vo, hon);
        kernel.set_out(v, hn);
        kernel.set_last(f, hn);
        
        HalfEdgeID hx = hon;
        
        GEL_ASSERT_EQ(kernel.face(hx), fo);
        while(kernel.face(hx) != f){
            kernel.set_face(hx, f);
            hx = kernel.next(hx);
        }
        
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
    
    VertexID Manifold::slit_vertex([[maybe_unused]] const VertexID vid, const HalfEdgeID h_in, const HalfEdgeID h_out)
    {
        return slit_one_ring(h_in, h_out);
    }

    VertexID Manifold::slit_one_ring(const HalfEdgeID h_in, const HalfEdgeID h_out)
    {
        // Check validity of input. Neither h_in nor h_out may be adjacent to a hole
        // i.e. their faces must be defined. They must also not be each others opposite
        // edges.
        if (kernel.face(h_in) == InvalidFaceID ||
            kernel.face(h_out) == InvalidFaceID ||
            kernel.opp(h_out) == h_in ||
            kernel.vert(kernel.opp(h_out)) != kernel.vert(h_in)
            )
            return InvalidVertexID;

        // Only one of h_in and h_out can be the opposite of a border.
        // Otherwise, there is no way of performing a manifold slit without immediately filling.
        if (kernel.next(kernel.opp(h_out)) == kernel.opp(h_in) &&
            kernel.face(kernel.opp(h_in)) == InvalidFaceID &&
            kernel.face(kernel.opp(h_out)) == InvalidFaceID)
            return InvalidVertexID; // can be replaced with an assert

        return slit_vertex_impl(h_in, h_out);
    }

    VertexID Manifold::slit_vertex_impl(const HalfEdgeID h_in, const HalfEdgeID h_out)
    {
        GEL_ASSERT_EQ(kernel.vert(kernel.opp(h_out)), kernel.vert(h_in));
        GEL_ASSERT_NEQ(h_in, kernel.opp(h_out));
        const VertexID v = kernel.vert(h_in);
        VertexID v_new = kernel.add_vertex();
        pos(v_new) = Vec3d(pos(v));

        // Go counter clockwise from h_out to h_in. Set all
        // halfedges that used to point to v to point to v_new
        HalfEdgeID h = kernel.prev(h_out);
        kernel.set_vert(h, v_new);
        while ( h != h_in) {
            h = kernel.prev(kernel.opp(h));
            kernel.set_vert(h, v_new);
        }

        HalfEdgeID h_in_opp = kernel.opp(h_in);
        HalfEdgeID h_new_in, h_new_in_opp;

        HalfEdgeID h_out_opp = kernel.opp(h_out);
        HalfEdgeID h_new_out,h_new_out_opp;

        // boundary case
        if (kernel.face(h_in_opp) == InvalidFaceID && kernel.face(h_out_opp) != InvalidFaceID)
        {
            h_new_in_opp = h_in_opp;
            h_new_in = kernel.prev(h_new_in_opp);
            h_in_opp = kernel.opp(h_new_in);
        } else {
            // Create two boundary edges that form a wedge between
            // h_in and h_in_opp -- if h_in_opp is not boundary
            h_new_in = kernel.add_halfedge();
            kernel.set_face(h_new_in, InvalidFaceID);
            glue(h_in_opp, h_new_in);

            h_new_in_opp = kernel.add_halfedge();
            kernel.set_face(h_new_in_opp, InvalidFaceID);
            glue(h_in, h_new_in_opp);

            link(h_new_in_opp, h_new_in);

            VertexID v_i = kernel.vert(h_in_opp);
            kernel.set_vert(h_new_in_opp, v_i);
            kernel.set_out(v_i, h_new_in);
        }

        // boundary case
        if (kernel.face(h_out_opp) == InvalidFaceID && kernel.face(h_in_opp) != InvalidFaceID) {
            h_new_out_opp = h_out_opp;
            h_new_out = kernel.next(h_new_out_opp);
            h_out_opp = kernel.opp(h_new_out);
        } else {
            // Create two boundary edges that for a wedge between
            // h_out and h_out_opp -- if h_out_opp is not boundary
            h_new_out = kernel.add_halfedge();
            kernel.set_face(h_new_out, InvalidFaceID);
            glue(h_out_opp, h_new_out);

            h_new_out_opp = kernel.add_halfedge();
            kernel.set_face(h_new_out_opp, InvalidFaceID);
            glue(h_out, h_new_out_opp);

            link(h_new_out, h_new_out_opp);

            VertexID v_o = kernel.vert(h_out);
            kernel.set_vert(h_new_out, v_o);
            kernel.set_out(v_o, h_new_out_opp);
        }

        link(h_new_out_opp, h_new_in_opp);
        link(h_new_in, h_new_out);

        kernel.set_vert(h_new_in, v);
        kernel.set_vert(h_new_out_opp, v_new);

        kernel.set_out(v, h_new_out);
        kernel.set_out(v_new, h_new_in_opp);

        return v_new;
    }

    
    HalfEdgeID Manifold::slit_edges(const HMesh::VertexSet& insel)
    {
        HalfEdgeID h;
        for(auto vid : insel)
        {
            HalfEdgeID h_in = InvalidHalfEdgeID, h_out = InvalidHalfEdgeID;
            Walker w = walker(vid);
            while(!w.full_circle())
            {
                if(insel.contains(w.vertex())) {
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
                VertexID v_new = slit_one_ring(h_in, h_out);
                h = walker(v_new).halfedge();
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
    }
    
    
    /**********************************************
     * Private functions
     **********************************************/
    

    
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
    
    void Manifold::merge(const Manifold& mergee) {
        Util::HashMap<FaceID,FaceID> fmap;
        Util::HashMap<HalfEdgeID,HalfEdgeID> hmap;
        Util::HashMap<VertexID,VertexID> vmap;

        fmap.reserve(kernel.allocated_faces());
        hmap.reserve(kernel.allocated_halfedges());
        vmap.reserve(kernel.allocated_vertices());
        
        for(auto v: mergee.vertices())
            vmap[v] = kernel.add_vertex();
        for(auto h: mergee.halfedges())
            hmap[h] = kernel.add_halfedge();
        for(auto f: mergee.faces())
            fmap[f] = kernel.add_face();
 
        for(auto f: mergee.faces()) {
            auto f_new = fmap[f];
            kernel.set_last(f_new, hmap[mergee.kernel.last(f)]);
        }
        
        for(auto h: mergee.halfedges()) {
            auto h_new = hmap[h];
            kernel.set_opp(h_new, hmap[mergee.kernel.opp(h)]);
            kernel.set_next(h_new, hmap[mergee.kernel.next(h)]);
            kernel.set_prev(h_new, hmap[mergee.kernel.prev(h)]);
            kernel.set_vert(h_new, vmap[mergee.kernel.vert(h)]);
            FaceID f = mergee.kernel.face(h);
            if (f == InvalidFaceID)
                kernel.set_face(h_new, InvalidFaceID);
            else
                kernel.set_face(h_new, fmap[f]);
        }
        
        for(auto v: mergee.vertices()) {
            auto v_new = vmap[v];
            pos(v_new) = mergee.pos(v);
            kernel.set_out(v_new, hmap[mergee.kernel.out(v)]);
        }
     }

    void Manifold::remove_face_if_degenerate(HalfEdgeID h)
    {
        // face is degenerate if there are only two halfedges in face loop
        if(kernel.next(kernel.next(h)) == h)
        {
            HalfEdgeID hn = kernel.next(h);
            HalfEdgeID ho = kernel.opp(h);
            HalfEdgeID hno = kernel.opp(hn);
            VertexID hv = kernel.vert(h);
            VertexID hnv = kernel.vert(hn);
            FaceID f = kernel.face(h);
            
            GEL_ASSERT_NEQ(ho, hn);
            GEL_ASSERT_NEQ(h , hno);
            GEL_ASSERT_NEQ(hv, hnv);
            
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

    VertexAttributeVector<size_t> build_manifold(
        Manifold& manifold,
        const std::span<Vec3d const> vertices, // need random access
        const std::span<size_t const> face_indexes, // need linear access
        const std::span<int const> face_vert_count
        )
    {
        // perform bound checking in advance for bad inputs
        const auto bound = vertices.size();
        for (const auto index: face_indexes) {
            GEL_ASSERT(index < bound, "build_manifold called with out of bounds input");
        }

        VertexAttributeVector<size_t> cluster_id;
        auto face_indexes_view = std::views::all(face_indexes);

        size_t running_total = 0;
        for (const int this_face_verts : face_vert_count) {
            // we need to break early if face_vert_count is too long
            if (running_total >= face_indexes.size())
                break;

            auto face_iter = face_indexes_view
                | std::views::take(this_face_verts)
                | std::views::transform([&](auto idx) { return vertices[idx]; });

            const FaceID f = manifold.add_face(face_iter);

            circulate_face_ccw(manifold, f, [&](VertexID v){
                cluster_id[v] = face_indexes_view.front();
                face_indexes_view = face_indexes_view | std::views::drop(1);
            });

            running_total += this_face_verts;
        }

        [[maybe_unused]] auto failed = stitch_mesh(manifold, cluster_id);
        // std::cerr << "unable to stitch: " << failed << "\n";
        IDRemap remap;
        manifold.cleanup(remap);
        cluster_id.cleanup(remap.vmap);
        return cluster_id;
    }

    VertexAttributeVector<size_t> build_manifold(
        Manifold& manifold,
        const std::span<Vec3d const> vertices, // need random access
        const std::span<size_t const> face_indexes, // need linear access
        const int face_vert_count // FIXME: C++23 adds ranges::repeat_view which means this overload can be removed with a template
        )
    {
        // perform bound checking in advance for bad inputs
        const auto bound = vertices.size();
        for (const auto index: face_indexes) {
            GEL_ASSERT(index < bound, "build_manifold called with out of bounds input");
        }

        VertexAttributeVector<size_t> cluster_id;
        auto face_indexes_view = std::views::all(face_indexes);

        for (size_t running_total = 0; running_total < face_indexes.size(); running_total += face_vert_count) {
            auto face_iter = face_indexes_view
                | std::views::take(face_vert_count)
                | std::views::transform([&](auto idx) { return vertices[idx]; });

            const FaceID f = manifold.add_face(face_iter);

            circulate_face_ccw(manifold, f, [&](VertexID v){
                cluster_id[v] = face_indexes_view.front();
                face_indexes_view = face_indexes_view | std::views::drop(1);
            });
        }

        [[maybe_unused]] auto failed = stitch_mesh(manifold, cluster_id);
        // std::cerr << "unable to stitch: " << failed << "\n";
        IDRemap remap;
        manifold.cleanup(remap);
        cluster_id.cleanup(remap.vmap);
        return cluster_id;
    }

    template<typename size_type, typename float_type, typename int_type>
    VertexAttributeVector<int_type> build_template(Manifold& m, size_type no_vertices,
                                  const float_type* vertvec,
                                  size_type no_faces,
                                  const int_type* facevec,
                                  const int_type* indices)
    {
        int k=0;
        VertexAttributeVector<int> cluster_id;
        for(int i=0;i<no_faces;++i) {
            vector<Vec3d> pts(facevec[i]);
            for(int j=0;j<facevec[i]; ++j) {
                const float_type* v = &vertvec[3*indices[j+k]];
                pts[j] = Vec3d(v[0], v[1], v[2]);
            }
            FaceID f = m.add_face(pts);
            int j=0;
            circulate_face_ccw(m, f, [&](VertexID v){
                cluster_id[v] = indices[j+k];
                ++j;
            });
            k += facevec[i];
        }
        stitch_mesh(m, cluster_id);
        IDRemap remap;
        m.cleanup(remap);
        cluster_id.cleanup(remap.vmap);
        return cluster_id;
    }
    
    VertexAttributeVector<int> build(Manifold& m, const TriMesh& mesh)
    {
        // A vector of 3's - used to tell build how many indices each face consists of
        vector<int> faces(mesh.geometry.no_faces(), 3);
        
        return build_template(m, static_cast<size_t>(mesh.geometry.no_vertices()),
                       reinterpret_cast<const float*>(&mesh.geometry.vertex(0)),
                       static_cast<size_t>(faces.size()),
                       static_cast<const int*>(&faces[0]),
                       reinterpret_cast<const int*>(&mesh.geometry.face(0)));
    }
    
    VertexAttributeVector<int> build(Manifold& m, size_t no_vertices,
                         const float* vertvec,
                         size_t no_faces,
                         const int* facevec,
                         const int* indices)
    {
        return build_template(m, no_vertices, vertvec, no_faces, facevec, indices);
    }
    
    VertexAttributeVector<int> build(Manifold& m, size_t no_vertices,
                         const double* vertvec,
                         size_t no_faces,
                         const int* facevec,
                         const int* indices)
    {
        return build_template(m, no_vertices, vertvec, no_faces, facevec, indices);
    }

    
    bool find_invalid_entities(const Manifold& m, VertexSet& vs, HalfEdgeSet& hs, FaceSet& fs)
    {
        bool valid = true;
        vs.clear();
        hs.clear();
        fs.clear();
        
        // Verify components of halfedges
        for(HalfEdgeID h : m.halfedges()){
            Walker j = m.walker(h);
            
            if(j.vertex() == InvalidVertexID){
                cout << "Halfedge lacks vert" << endl;
                hs.insert(h);
                valid = false;
            }
            if(j.next().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks next" << endl;
                hs.insert(h);
                valid = false;
            }
            if(j.prev().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks prev" << endl;
                hs.insert(h);
                valid = false;
            }
            if(j.opp().halfedge() == InvalidHalfEdgeID){
                cout << "Halfedge lacks opp" << endl;
                hs.insert(h);
                valid = false;
            }
            
        }
        // Verify components of vertices
        for(VertexID v : m.vertices()){
            vector<VertexID> link;
            
            // circulate the halfedges of vertex
            for(Walker j = m.walker(v); !j.full_circle(); j = j.circulate_vertex_cw()){
                // test halfedges around v
                if(j.halfedge() == InvalidHalfEdgeID){
                    cout << "Vertex circulation produced invalid halfedge" << endl;
                    valid = false;
                    vs.insert(v);
                    break;
                }
                VertexID ring_v = j.vertex();
                if(!m.in_use(ring_v))
                {
                    cout << "Invalid vertex: " << ring_v << " in one-ring of vertex" << endl;
                    valid = false;
                    vs.insert(v);
                    break;
                }
                
                // test one-ring for multiple occurences of vertex
                if(find(link.begin(), link.end(), ring_v) != link.end()){
                    cout << "Vertex appears two times in one-ring of vertex" << endl;
                    valid = false;
                    vs.insert(v);
                    break;
                }
                link.push_back(ring_v);
                
                // test for infinite loop around vertex
                if(static_cast<size_t>(j.no_steps()) > m.no_vertices()){
                    cout << "Vertex loop CW contains more vertices than manifold" << endl;
                    valid = false;
                    vs.insert(v);
                    break;
                }
            }
            
            for(Walker j = m.walker(v); !j.full_circle(); j = j.circulate_vertex_ccw()) {
                if(static_cast<size_t>(j.no_steps()) > m.no_vertices()){
                    cout << "Vertex loop CCW contains more vertices than manifold" << endl;
                    valid = false;
                    vs.insert(v);
                    break;
                }
            }
            
            if(link.size()==1) {
                cout << "Vertex contains only a single incident edge" << endl;
                vs.insert(v);
                valid = false;
            }
        }
        // verify components of faces
        for(FaceID f : m.faces()){
            // count edges on face
            Walker j = m.walker(f);
            
            vector<VertexID> face_verts;
            for(; !j.full_circle(); j = j.circulate_face_cw()){
                if(find(begin(face_verts), end(face_verts), j.vertex()) == face_verts.end())
                    face_verts.push_back(j.vertex());
                else {
                    cout << "Face contains same vertex twice" << endl;
                    valid = false;
                    fs.insert(f);
                    break;
                }
                     
                // test that all halfedges in faces bind properly to their face
                if(j.face() != f){
                    cout << "Face is inconsistent, halfedge is not bound to face" << endl;
                    valid = false;
                    fs.insert(f);
                    break;
                }
            }
            // test faces for valid geometrical properties
            if(j.no_steps() < 3){
                cout << "Face contains less than 3 edges" << endl;
                fs.insert(f);
                valid = false;
            }
            // test for infinite loop around face
            if(j.no_steps() > m.no_halfedges() * 0.5f){
                cout << "Face loop contains more halfedges than manifold" << endl;
                fs.insert(f);
                valid = false;
            }
        }
        return valid;
    }
    
    bool valid(const Manifold& m)
    {
        VertexSet vs;
        HalfEdgeSet hs;
        FaceSet fs;
        return find_invalid_entities(m, vs, hs, fs);
    }
    
    void bbox(const Manifold& m, Manifold::Vec& pmin, Manifold::Vec& pmax) {
        if(m.no_vertices()==0)
            return;
        
        pmin = Vec3d(DBL_MAX);
        pmax = Vec3d(-DBL_MAX);
        
        for (auto v: m.vertices()) {
            pmin = v_min(m.pos(v), pmin);
            pmax = v_max(m.pos(v), pmax);
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
    
    
    
    bool Manifold::precond_collapse_edge(HalfEdgeID h) const
    {
        Walker hew = walker(h);
        HalfEdgeID ho = hew.opp().halfedge();
        VertexID v0 = hew.opp().vertex();
        VertexID v1 = hew.vertex();
        
        // get the one-ring of v0
        vector<VertexID> link0;
        vector<FaceID> faces0;
        int k = 0;
        for(Walker vj = walker(h);
            !vj.full_circle(); vj = vj.circulate_vertex_ccw()){
            link0.push_back(vj.vertex());
            if(vj.halfedge() != h)
                faces0.push_back(vj.face());
            if(++k>no_vertices())
            {
                cout << "precond_collapse failed: mesh is corrupted" << endl;
                return false;
            }
        }
        if(link0.size() <= 1) {
            cout << "precond_collapse failed: link0 size is <=1" << endl;
            return false;
        }
        
        // get the one-ring of v1
        vector<VertexID> link1;
        vector<FaceID> faces1;
        k=0;
        for(Walker vj = walker(ho);
            !vj.full_circle(); vj = vj.circulate_vertex_ccw()){
            link1.push_back(vj.vertex());
            if(vj.halfedge() != ho)
                faces1.push_back(vj.face());
            if(++k>no_vertices())
            {
                cout << "precond_collapse failed: mesh is corrupted" << endl;
                return false;
            }
        }
        if(link1.size() <= 1) {
            cout << "precond_collapse failed: link1 size is <=1" << endl;
            return false;
        }

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
        if(fisect.size() > 0) {
//            cout << "precond_collapse failed: same face " << fisect[0] << " in both 1-rings" << endl;
            return false;
        }
        
         k = 0;
        // if the adjacent face is a triangle (see 2)
        if(hew.next().next().next().halfedge() == h){
            VertexID v = hew.next().vertex();
            
            // valency test (see 5)
            if(valency(v) < 3)
            {
//                cout << "precond_collapse failed: left vertex in triangle has val<3" << endl;
                return false;
            }

            
            // remove the vertex shared by the two rings from the intersection
            vector<VertexID>::iterator iter;
            iter = find(lisect.begin(), lisect.end(), v);
            if (iter != lisect.end()) {
                lisect.erase(iter);
                ++k;
            }
        }
        // if the adjacent face is a triangle (see 2)
        if(hew.opp().next().next().next().halfedge() == hew.opp().halfedge()){
            VertexID v = hew.opp().next().vertex();
            
            // valency test (see 5)
            if(valency(v) < 3)
            {
//                cout << "precond_collapse failed: right vertex in triangle has val<3" << endl;
                return false;
            }

            
            // remove the vertex shared by the two rings from the intersection
            vector<VertexID>::iterator iter;
            iter = find(lisect.begin(), lisect.end(), v);
            if (iter != lisect.end()) {
                lisect.erase(iter);
                ++k;
            }
        }
        // double edge test (see 3)
        if(lisect.size() != 0)
        {
//            cout << "precond_collapse failed: vertex shared by both 1-rings" << endl;
            return false;
        }
        
//        // tetrahedon test (see 4)
//        if(k == 2 && (link0.size() + link1.size() == 6))
//        {
//            cout << "precond_collapse failed: tet test" << endl;
//            return false;
//        }
        
        // test that we do not merge holes (see 6)
        if(boundary(v0) && boundary(v1) && !boundary(h))
        {
//            cout << "precond_collapse failed: would merge holes" << endl;
            return false;
        }
        
        return true;
    }

    
    bool Manifold::precond_flip_edge(HalfEdgeID h) const 
    {
        Walker j = walker(h);
        
        FaceID hf = j.face();
        FaceID hof = j.opp().face();
        
        // boundary case
        if(hf == InvalidFaceID || hof == InvalidFaceID)
           return false;
        
        
        // We can only flip an edge if both incident polygons are triangles.
        if(no_edges(hf) != 3 || no_edges(hof) !=3)
            return false;
        

        // non boundary vertices with a valency of less than 4(less than 3 after operation) degenerates mesh.
        VertexID hv = j.vertex();
        VertexID hov = j.opp().vertex();
        if((valency(hv) < 4 && !boundary(hv)) || (valency(hov) < 4 && !boundary(hov))){
            return false;
        }
        
        // Disallow flip if vertices being connected already are.
        VertexID hnv = j.next().vertex();
        VertexID honv = j.opp().next().vertex();
        if(connected(hnv, honv)){
           return false;
        }
        
        return true;
    }

    Manifold::Vec Manifold::normal(VertexID v) const
    {
        Manifold::Vec p0 = positions[v];
        Manifold::Vec n(0);
        circulate_vertex_ccw(*this, v, [&](Walker& w) {
            if(w.face() != InvalidFaceID) {
                Manifold::Vec e0 = cond_normalize(positions[w.prev().opp().vertex()] - p0);
                Manifold::Vec e1 = cond_normalize(positions[w.vertex()] - p0);
                Manifold::Vec n_face = cond_normalize(cross(e1, e0));
                n += acos(max(-1.0, min(1.0, dot(e0, e1)))) * n_face;
            }
        });
        return cond_normalize(n);
    }
    
    Manifold::Vec Manifold::area_normal(FaceID f) const
    {
        GEL_ASSERT_NEQ(f, InvalidFaceID);
        std::vector<Manifold::Vec> vertices;
        Vec c(0.0);
        int k = circulate_face_ccw(*this, f, [&](VertexID vid) {
            Vec p = positions[vid];
            c += p;
            vertices.push_back(p);
        });
        c /= k;
        Manifold::Vec norm(0);
        for(int i=0;i<k;++i)
            norm += cross(vertices[i]-c,vertices[(i+1)%k]-c);
        return 0.5 * norm;
    }
    
    double Manifold::area(FaceID fid) const
    {
        // Get all projected vertices
        std::vector<Manifold::Vec> vertices;
        int N = circulate_face_ccw(*this, fid, [&](VertexID vid) {
            vertices.push_back(positions[vid]);
        });
        
        double area = 0;
        Manifold::Vec norm = normal(fid);
        for(int i = 1; i < N-1; ++i)
            area += 0.5 * dot(norm,cross(vertices[i]-vertices[0], vertices[(i+1 )]-vertices[0]));
        return area;
    }


    bool closed(const Manifold& m)
    {
        for(auto h: m.halfedges())
            if(m.walker(h).face() == InvalidFaceID)
                return false;
        return true;
    }

    Manifold::Vec barycenter(const Manifold& m) {
        Manifold::Vec c(0);
        for (auto v: m.vertices()) {
            c += m.pos(v);
        }
        return c/m.no_vertices();
    }

    double area(const Manifold& m) {
        double a = 0;
        for (auto f: m.faces())
            a += area(m, f);
        return a;
    }

    double volume(const Manifold& _m) {
        Manifold m = Manifold(_m);
        triangulate(m);
        Vec3d c;
        float r;
        bsphere(m, c, r);
        double vol = 0;
        for (auto f: m.faces()) {
            Vec3d p[3];
            int i=0;
            for(auto v: m.incident_vertices(f))
                p[i++] = m.pos(v);
            vol += determinant(Mat3x3d(p[0]-c, p[1]-c, p[2]-c));
        }
        return vol / 6.0;
    }


    vector<HMesh::Manifold> connected_components(const HMesh::Manifold& m) {
        auto visited = FaceAttributeVector<int>(m.allocated_faces(), 0);
        vector<HMesh::Manifold> components;
        for (auto face: m.faces()) {
            if (visited[face] == 0) {
                stack<FaceID> face_stack;
                face_stack.push(face);
                visited[face] = 1;
                components.push_back(Manifold());
                Manifold& comp_mesh = components.back();
                VertexAttributeVector<int> orig_id;
                while(!face_stack.empty()) {
                    FaceID current_face = face_stack.top();
                    face_stack.pop();
                    vector<Vec3d> pts;
                    vector<int> ids;
                    for (auto v: m.incident_vertices(current_face)) {
                        pts.push_back(m.pos(v));
                        ids.push_back(v.get_index());
                        for (auto neighbor: m.incident_faces(v))
                            if (m.in_use(neighbor) and visited[neighbor] == 0) {
                                face_stack.push(neighbor);
                                visited[neighbor] = 1;
                            }
                    }
                    FaceID f_new = comp_mesh.add_face(pts);
                    int i=0;
                    for (auto v: comp_mesh.incident_vertices(f_new))
                        orig_id[v] = ids[i++];
                }
                stitch_mesh(comp_mesh, orig_id);
                comp_mesh.cleanup();
            }
        }
        return components;
    }


    int count_boundary_curves(const HMesh::Manifold& m) {
        int num = 0;
        HalfEdgeAttributeVector<int> visited(m.allocated_halfedges(), 0);
        for (auto h0: m.halfedges()) {
            if (visited[h0] == 1)
                continue;
            Walker w = m.walker(h0);
            if (!m.in_use(w.face())) {
                num += 1;
                while (visited[w.halfedge()] == 0) {
                    visited[w.halfedge()] = 1;
                    w = w.next();
                }
            }
        }
        return num;
    }


}
