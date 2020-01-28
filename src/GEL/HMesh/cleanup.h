/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file caps_and_needles.h
 * @brief Simple tools for improving polygonal meshes by removing bad triangles.
 */

#ifndef __HMESH_CAPS_AND_NEEDLES_H__
#define __HMESH_CAPS_AND_NEEDLES_H__

namespace HMesh
{
    class Manifold;
    template<class T> class VertexAttributeVector;
    /** \brief Remove caps from a manifold consisting of only triangles.
    A cap is a triangle with two very small angles and an angle close to pi, however a cap does not necessarily have a very short edge.
    Set the ang_thresh to a value close to pi. The closer to pi the _less_ sensitive the cap removal.
    A cap is removed by flipping the (long) edge E opposite to the vertex V with the angle close to pi. 
    However, the function is more complex. Read code and document more carefully !!! */
    void remove_caps(Manifold& m, float thresh);

    /** \brief Remove needles from a manifold consisting of only triangles.
    A needle is a triangle with a single very short edge. It is moved by collapsing the short edge. 
    The thresh parameter sets the length threshold as a fraction of the average edge length.		 */
    void remove_needles(Manifold& m, float thresh=0.1, bool averagePositions = true);
    
    /** \brief Stitch together edges whose endpoints coincide geometrically. 
     This function allows you to create a mesh as a bunch of faces and then stitch these together
     to form a coherent whole. What this function adds is a spatial data structure to find out
     which vertices coincide. The return value is the number of edges that could not be stitched. 
     Often this is because it would introduce a non-manifold situation.*/
    int stitch_mesh(Manifold& m, double rad);
    int stitch_mesh(Manifold& m, const VertexAttributeVector<int>& cluster_id);

    /** \brief Stitches the mesh together, splits edges that could not be stitched and goes again.
     This function thereby handles situations where stitch mesh would not have worked. */
    void stitch_more(Manifold& m, double rad);

    /** \brief This function replaces holes by faces.
     It is really a simple function that just finds all loops of edges next to missing faces.
     You can specify, max_size, the maximum hole size to close. */
    void close_holes(Manifold& m, int max_size=100);
    
    /** \brief Flip the orientation of a mesh.
     After calling this function, normals will point the other way and clockwise becomes 
     counter clockwise */
    void flip_orientation(Manifold& m);
    
    /** Remove valence two vertices. */
    void remove_valence_two_vertices(Manifold & m);

    /** This function merges pairs of boundary vertices, provided there are exactly two such vertices
     at a given point in space, that they are not in each others' one ring and  that the one rings are disjoint. */
    void merge_coincident_boundary_vertices(Manifold& m, double rad=1e-30);
}

#endif
