/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
@file Walker.h
 Contains class for walking a mesh.
 */

#pragma once

#include "Manifold.h"

namespace HMesh
{
    /** Class for traversing the entities of a HMesh. This class can work as 
     both a vertex and a face circulator but also move in other ways. It is,
     in fact, the only way to traverse the mesh from the users point of view. */
    class Walker
    {
    public:
        /// construct from kernel and a halfedge
        Walker(const ConnectivityKernel& _ck, HalfEdgeID _current);

        /// returned walker has made one step to the next halfedge 
        Walker next() const;
        /// returned walker has made one step to the previous halfedge 
        Walker prev() const;
        /// returned walker has made one step to the opposite halfedge 
        Walker opp() const;

        /// returned walker has circulated vertex clockwise one step
        Walker circulate_vertex_cw() const;
        /// returned walker has circulated vertex counterclockwise one step
        Walker circulate_vertex_ccw() const;

        /// returned walker has circulated face clockwise one step
        Walker circulate_face_cw() const;
        /// returned walker has circulated face counterclockwise one step
        Walker circulate_face_ccw() const;

        /// test if the walker has reached its initial halfedge
        bool full_circle() const;
        
        /// number of steps taken
        int no_steps() const;

        /// get ID of vertex pointed to by current halfedge of walker
        VertexID vertex() const; 
        /// get ID of face owning current halfedge of walker
        FaceID face() const; 
        /// get ID of current halfedge of walker
        HalfEdgeID halfedge() const;
        
        /// assignment operator
        Walker operator =(const Walker& w);

    private:
        const ConnectivityKernel* ck;
        HalfEdgeID current;
        HalfEdgeID last;
        int counter;

        Walker(const ConnectivityKernel& _ck, HalfEdgeID _current, HalfEdgeID _last, int _counter);
    };

    inline Walker::Walker(const ConnectivityKernel& _ck, HalfEdgeID _current) 
        : ck(&_ck), current(_current), last(_current), counter(0){}

    inline Walker::Walker(const ConnectivityKernel& _ck, HalfEdgeID _current, HalfEdgeID _last, int _counter)
        : ck(&_ck), current(_current), last(_last), counter(_counter){}

    inline Walker Walker::next() const
    { return Walker(*ck, ck->next(current), last, counter + 1); }

    inline Walker Walker::prev() const
    { return Walker(*ck, ck->prev(current), last, counter + 1); }

    inline Walker Walker::opp() const
    { return Walker(*ck, ck->opp(current), last, counter + 1); }

    inline Walker Walker::circulate_vertex_cw() const
    { return Walker(*ck, ck->next(ck->opp(current)), last, counter + 1); }

    inline Walker Walker::circulate_vertex_ccw() const
    { return Walker(*ck, ck->opp(ck->prev(current)), last, counter + 1); }

    inline Walker Walker::circulate_face_cw() const
    { return Walker(*ck, ck->prev(current), last, counter + 1); }

    inline Walker Walker::circulate_face_ccw() const
    { return Walker(*ck, ck->next(current), last, counter + 1); }

    inline bool Walker::full_circle() const
    { return (counter > 0 && current == last) ? true : false; }
	
    inline int Walker::no_steps() const
    { return counter; }

    inline VertexID Walker::vertex() const
    { return ck->vert(current); }

    inline FaceID Walker::face() const
    { return ck->face(current); }

    inline HalfEdgeID Walker::halfedge() const
    { return current; }

    inline Walker Walker::operator =(const Walker& w)
    { 
        current = w.current;
        counter = w.counter;
        ck = w.ck;
        last = w.last;
        return *this;
    }
    
 

}

