/**
 * @file triangulate.h
 * @brief Triangulating the faces of a mesh.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_TRIANGULATE__H
#define __HMESH_TRIANGULATE__H

#include "Manifold.h"

namespace HMesh
{
    /// Naive division of polygons into triangles.
    void triangulate_by_edge_face_split(Manifold& m);

    /// Try to respect curvature to create a better triangulation.
    void curvature_triangulate(Manifold& m);

    /// Naive triangulation by connecting to center point.
    void triangulate_by_vertex_face_split(Manifold& m);

    /// Triangulate by connecting the points forming the shortest edge.
    void shortest_edge_triangulate(Manifold& m);

    /** \brief Triangulate a polygonal face by repeatedly calling split_face.
    split_face_triangulate iteratively splits triangles off a polygon. 
    The first triangle split off is the one connecting f.last().vert() and f.last().next().next().vert(). */
    void triangulate_face_by_edge_split(Manifold& m, FaceID f);


}

#endif