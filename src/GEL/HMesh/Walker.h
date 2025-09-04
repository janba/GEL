/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */
#pragma once
#ifndef GEL_HMESH_WALKER_H
#define GEL_HMESH_WALKER_H

/**
@file Walker.h
 Contains class for walking a mesh.
 */


#include <functional>

namespace HMesh
{
/// Class for traversing the entities of a HMesh. This class can work as
/// both a vertex and a face circulator but also move in other ways. It is,
/// in fact, the only way to traverse the mesh from the user's point of view.
class Walker {
public:
    /// construct from kernel and a halfedge
    Walker(const ConnectivityKernel& _ck, HalfEdgeID _current);
    // Compiler requires a copy constructor when Walker is passed by value
    Walker(const Walker& other) = default;

    /// @return new walker that made one step to the next halfedge
    [[nodiscard]] Walker next() const;
    /// @return new walker that made one step to the previous halfedge
    [[nodiscard]] Walker prev() const;
    /// @return new walker that made one step to the opposite halfedge
    [[nodiscard]] Walker opp() const;

    /// @return new walker that circulated vertex clockwise one step
    [[nodiscard]] Walker circulate_vertex_cw() const;
    /// @return new walker that circulated vertex counterclockwise one step
    [[nodiscard]] Walker circulate_vertex_ccw() const;

    /// @return new walker that circulated face clockwise one step
    [[nodiscard]] Walker circulate_face_cw() const;
    /// @return new walker that circulated face counterclockwise one step
    [[nodiscard]] Walker circulate_face_ccw() const;

    /// test if the walker has reached its initial halfedge
    [[nodiscard]] bool full_circle() const;

    /// number of steps taken
    [[nodiscard]] int no_steps() const;

    /// get ID of vertex pointed to by the current halfedge of the walker
    [[nodiscard]] VertexID vertex() const;
    /// get ID of face owning current halfedge of walker
    [[nodiscard]] FaceID face() const;
    /// get ID of current halfedge of walker
    [[nodiscard]] HalfEdgeID halfedge() const;
    /// Get ID of either halfedge or ID - whichever has the smaller index.
    [[nodiscard]] HalfEdgeID hmin() const;

    Walker& operator =(const Walker& w);

private:
    std::reference_wrapper<const ConnectivityKernel> ck;
    HalfEdgeID current;
    HalfEdgeID last;
    int counter;

    Walker(const ConnectivityKernel& _ck, HalfEdgeID _current, HalfEdgeID _last, int _counter);
};

inline Walker::Walker(const ConnectivityKernel& _ck, const HalfEdgeID _current)
    : ck(_ck), current(_current), last(_current), counter(0) {}

inline Walker::Walker(const ConnectivityKernel& _ck, const HalfEdgeID _current, const HalfEdgeID _last,
                      const int _counter)
    : ck(_ck), current(_current), last(_last), counter(_counter) {}

inline Walker Walker::next() const
{
    return {ck, ck.get().next(current), last, counter + 1};
}

inline Walker Walker::prev() const
{
    return {ck, ck.get().prev(current), last, counter + 1};
}

inline Walker Walker::opp() const
{
    return {ck, ck.get().opp(current), last, counter + 1};
}

inline Walker Walker::circulate_vertex_cw() const
{
    return {ck, ck.get().next(ck.get().opp(current)), last, counter + 1};
}

inline Walker Walker::circulate_vertex_ccw() const
{
    return {ck, ck.get().opp(ck.get().prev(current)), last, counter + 1};
}

inline Walker Walker::circulate_face_cw() const
{
    return {ck, ck.get().prev(current), last, counter + 1};
}

inline Walker Walker::circulate_face_ccw() const
{
    return {ck, ck.get().next(current), last, counter + 1};
}

inline bool Walker::full_circle() const
{
    return (counter > 0 && current == last);
}

inline int Walker::no_steps() const
{
    return counter;
}

inline VertexID Walker::vertex() const
{
    return ck.get().vert(current);
}

inline FaceID Walker::face() const
{
    return ck.get().face(current);
}

inline HalfEdgeID Walker::halfedge() const
{
    return current;
}


inline HalfEdgeID Walker::hmin() const
{
    return (current < ck.get().opp(current)) ? current : ck.get().opp(current);
}


inline Walker& Walker::operator=(const Walker& w) = default;
}

#endif
