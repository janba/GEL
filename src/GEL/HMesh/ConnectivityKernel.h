/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file ConnectivityKernel.h
 * @brief The data structure under HMesh.
 */

#ifndef __HMESH_CONNECTIVITY_KERNEL_H__
#define __HMESH_CONNECTIVITY_KERNEL_H__

#include <vector>
#include <map>
#include "ItemVector.h"
#include "ItemID.h"
#include "Iterators.h"

namespace HMesh
{
    struct Vertex;
    struct Face;
    struct HalfEdge;

    typedef ItemVector<Vertex>::IDType VertexID;
    typedef ItemVector<Face>::IDType FaceID;
    typedef ItemVector<HalfEdge>::IDType HalfEdgeID;

    /// The vertex struct. This contains just a single outgoing halfedge.
    struct Vertex
    { 
        HalfEdgeID out; 
    };

    /// The face struct. Contains a single halfedge
    struct Face
    { 
        HalfEdgeID last; 
    };

    /// The halfedge struct. Contains IDs of next, previous, and opposite edges as well as incident face and vertex.
    struct HalfEdge
    {
        HalfEdgeID next;
        HalfEdgeID prev;
        HalfEdgeID opp;
        VertexID vert;
        FaceID face;
    };
    
    
    typedef IDIterator<Vertex> VertexIDIterator;
    typedef IDIterator<Face> FaceIDIterator;
    typedef IDIterator<HalfEdge> HalfEdgeIDIterator;

    static const VertexID InvalidVertexID;
    static const FaceID InvalidFaceID;
    static const HalfEdgeID InvalidHalfEdgeID;
    
    typedef std::map<VertexID, VertexID> VertexIDRemap;
    typedef std::map<FaceID, FaceID> FaceIDRemap;
    typedef std::map<HalfEdgeID, HalfEdgeID> HalfEdgeIDRemap;
  
    /// The IDRemap struct is just used for garbage collection.
    struct IDRemap
    {
        VertexIDRemap vmap;
        FaceIDRemap fmap;
        HalfEdgeIDRemap hmap;
    };

    /** The connectivity kernel is basically an aggregate of ItemVectors for vertices, faces, and halfedges.
     This class contains no geometry information - only information about connectivitiy. Arguably it abstracts
     away the implementation from the ConnectivityKernel class making it possible, for instance, to use a different kernel. */
    class ConnectivityKernel
    {
    public:
        
        /// number of active vertices in kernel
        size_t no_vertices() const;
        /// number of active faces in kernel
        size_t no_faces() const;
        /// number of active halfedges in kernel
        size_t no_halfedges() const;
        
        /// number of total vertices in kernel
        size_t allocated_vertices() const;
        /// number of total faces in kernel
        size_t allocated_faces() const;
        /// number of total halfedges in kernel
        size_t allocated_halfedges() const;
        
        /// check if ID of vertex is in use
        bool in_use(VertexID id) const;
        /// check if ID of face is in use
        bool in_use(FaceID id) const;
        /// check if ID of halfedge is in use
        bool in_use(HalfEdgeID id) const;
        
        /// Iterator to first VertexID, optional argument defines if unused items should be skipped
        VertexIDIterator vertices_begin(bool skip = true) const;
        /// Iterator to first FaceID, optional argument defines if unused items should be skipped
        FaceIDIterator faces_begin(bool skip = true) const;
        /// Iterator to first HalfEdgeID, optional argument defines if unused items should be skipped
        HalfEdgeIDIterator halfedges_begin(bool skip = true) const;
        
        /// Iterator to past the end VertexID
        VertexIDIterator vertices_end() const;
        /// Iterator topast the end FaceID
        FaceIDIterator faces_end() const;
        /// Iterator to past the end HalfEdgeID
        HalfEdgeIDIterator halfedges_end() const;

        /// get the ID of next halfedge, given by current halfedge ID
        HalfEdgeID next(HalfEdgeID current) const;
        /// get the ID of previous halfedge, given by current halfedge ID
        HalfEdgeID prev(HalfEdgeID current) const;
        /// get the ID of opposite halfedge, given by current halfedge ID
        HalfEdgeID opp(HalfEdgeID current) const;
        /// get the ID of outgoing halfedge, given by current vertex ID
        HalfEdgeID out(VertexID current) const;
        /// get the ID of last halfedge of current face ID
        HalfEdgeID last(FaceID current) const;
        /// get the ID of vertex pointed to by current halfedge ID
        VertexID vert(HalfEdgeID id) const;
        /// get the ID of face owning current halfedge ID
        FaceID face(HalfEdgeID id) const;

        /// add vertex to kernel
        VertexID add_vertex();
        /// add face to kernel
        FaceID add_face();
        /// add halfedge to kernel
        HalfEdgeID add_halfedge();

        /// remove vertex from kernel, given by ID
        void remove_vertex(VertexID id);
        /// remove face from kernel, given by ID
        void remove_face(FaceID id);
        /// remove halfedge from kernel, given by ID
        void remove_halfedge(HalfEdgeID id);

        /// set the ID of next halfedge of current halfedge to next
        void set_next(HalfEdgeID current, HalfEdgeID next);
        /// set the ID of previous halfedge of current halfedge to prev
        void set_prev(HalfEdgeID current, HalfEdgeID prev);
        /// set the ID of opposite halfedge of current halfedge to opp
        void set_opp(HalfEdgeID current, HalfEdgeID opp);
        /// set the ID of outgoing halfedge of current vertex to out
        void set_out(VertexID current, HalfEdgeID out);
        /// set the ID of last halfedge of current face to last
        void set_last(FaceID current, HalfEdgeID last);
        /// set the ID of vertex pointed to by current halfedge to vert
        void set_vert(HalfEdgeID current, VertexID vert);
        /// set the ID of face owning current halfedge to face
        void set_face(HalfEdgeID current, FaceID face);

        /// Clean up unused space in vectors - WARNING! Invalidates existing handles!
        void cleanup(IDRemap& map);

        /// clear the kernel
        void clear();
        
    private:

        ItemVector<Vertex> vertices;
        ItemVector<Face> faces;
        ItemVector<HalfEdge> halfedges;
    };

    inline VertexID ConnectivityKernel::add_vertex()
    { return vertices.add(Vertex()); }

    inline FaceID ConnectivityKernel::add_face()
    { return faces.add(Face()); }

    inline HalfEdgeID ConnectivityKernel::add_halfedge()
    { return halfedges.add(HalfEdge()); }

    inline void ConnectivityKernel::remove_vertex(VertexID id)
    { vertices.remove(id); }

    inline void ConnectivityKernel::remove_face(FaceID id)
    { faces.remove(id); }

    inline void ConnectivityKernel::remove_halfedge(HalfEdgeID id)
    { halfedges.remove(id); }


    inline HalfEdgeID ConnectivityKernel::next(HalfEdgeID id) const
    { return halfedges[id].next; }

    inline HalfEdgeID ConnectivityKernel::prev(HalfEdgeID id) const
    { return halfedges[id].prev; }

    inline HalfEdgeID ConnectivityKernel::opp(HalfEdgeID id) const
    { return halfedges[id].opp; }

    inline HalfEdgeID ConnectivityKernel::out(VertexID id) const
    { return vertices[id].out; }

    inline HalfEdgeID ConnectivityKernel::last(FaceID id) const
    { return faces[id].last; }

    inline VertexID ConnectivityKernel::vert(HalfEdgeID id) const
    { return halfedges[id].vert; }

    inline FaceID ConnectivityKernel::face(HalfEdgeID id) const
    { return halfedges[id].face; }

    inline void ConnectivityKernel::set_next(HalfEdgeID id, HalfEdgeID next)
    { halfedges[id].next = next; }

    inline void ConnectivityKernel::set_prev(HalfEdgeID id, HalfEdgeID prev)
    { halfedges[id].prev = prev; }

    inline void ConnectivityKernel::set_opp(HalfEdgeID id, HalfEdgeID opp)
    {halfedges[id].opp = opp; }

    inline void ConnectivityKernel::set_out(VertexID id, HalfEdgeID out)
    { vertices[id].out = out; }

    inline void ConnectivityKernel::set_last(FaceID id, HalfEdgeID last)
    { faces[id].last = last; }

    inline void ConnectivityKernel::set_vert(HalfEdgeID id, VertexID vert)
    { halfedges[id].vert = vert; }

    inline void ConnectivityKernel::set_face(HalfEdgeID id, FaceID face)
    { halfedges[id].face = face; }



    inline size_t ConnectivityKernel::no_vertices() const
    { return vertices.size(); }

    inline size_t ConnectivityKernel::no_faces() const
    { return faces.size(); }

    inline size_t ConnectivityKernel::no_halfedges() const
    { return halfedges.size(); }



    inline size_t ConnectivityKernel::allocated_vertices() const
    { return vertices.allocated_size(); }

    inline size_t ConnectivityKernel::allocated_faces() const
    { return faces.allocated_size(); }

    inline size_t ConnectivityKernel::allocated_halfedges() const
    { return halfedges.allocated_size(); }



    inline bool ConnectivityKernel::in_use(VertexID id) const
    { return vertices.in_use(id); }

    inline bool ConnectivityKernel::in_use(FaceID id) const
    { return faces.in_use(id); }

    inline bool ConnectivityKernel::in_use(HalfEdgeID id) const
    { return halfedges.in_use(id); }

    inline VertexIDIterator ConnectivityKernel::vertices_begin(bool skip) const
    { return VertexIDIterator(vertices, vertices.index_begin(skip), skip); }
    inline  FaceIDIterator ConnectivityKernel::faces_begin(bool skip) const
    { return FaceIDIterator(faces, faces.index_begin(skip), skip); }
    inline HalfEdgeIDIterator ConnectivityKernel:: halfedges_begin(bool skip) const
    { return HalfEdgeIDIterator(halfedges, halfedges.index_begin(skip), skip); }
    
    
    inline VertexIDIterator ConnectivityKernel::vertices_end() const
    { return VertexIDIterator(vertices, vertices.index_end()); }
    inline FaceIDIterator ConnectivityKernel::faces_end() const
    { return FaceIDIterator(faces, faces.index_end()); }
    inline HalfEdgeIDIterator ConnectivityKernel::halfedges_end() const
    { return HalfEdgeIDIterator(halfedges, halfedges.index_end()); }
    
    
    inline void ConnectivityKernel::clear()
    {
        vertices.clear();
        faces.clear();
        halfedges.clear();
    }
}

#endif