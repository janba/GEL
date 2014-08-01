/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Manifold.h
 * @brief The Manifold class is the main data structure of HMesh - the actual mesh.
 */

#pragma once

#include <algorithm>
#include "../CGLA/Vec3d.h"

#include "ConnectivityKernel.h"
#include "Iterators.h"
#include "Walker.h"
#include "AttributeVector.h"


namespace Geometry
{
    // forward declaration
    class TriMesh;
    class IndexedFaceSet;
}

namespace HMesh
{
    /** The Manifold class represents a halfedge based mesh. Since meshes based on the halfedge
     representation must be manifold (although exceptions could be made) the class is thus named.
     Manifold contains many functions for mesh manipulation and associated the position attribute
     with vertices.*/
    
    class Manifold
    {
    public:
        
        /// Vector type used for positions of vertices.
        typedef CGLA::Vec3d Vec;
        
        /// Default constructor
        Manifold();

        /** \brief Build a manifold. 
        The arguments are the number of vertices, no_vertices, the vector of vertices, vertvec, the number of faces, no_faces. 
        facevec is an array where each entry indicates the number of vertices in that face.
        The array indices contains all the corresponding vertex indices in one concatenated list. */
        void build( size_t no_vertices,
                    const float* vertvec,
                    size_t no_faces,
                    const int* facevec,
                    const int* indices);

        /** \brief Build a manifold.
         This function is for vertices given in double precision.
         The arguments are the number of vertices, no_vertices, the vector of vertices, vertvec, the number of faces, no_faces.
         facevec is an array where each entry indicates the number of vertices in that face.
         The array indices contains all the corresponding vertex indices in one concatenated list. */
        void build( size_t no_vertices,
                   const double* vertvec,
                   size_t no_faces,
                   const int* facevec,
                   const int* indices);

        /// Build a manifold from a TriMesh
        void build(const Geometry::TriMesh& mesh);
        
        /** Add a face to the Manifold.
         This function is provided a vector of points in space and transforms it into a single 
         polygonal face calling build. It is purely for convenience. */
        FaceID add_face(std::vector<Manifold::Vec> points);

        /** Removes a face from the Manifold. If it is an interior face it is simply replaces
         by an InvalidFaceID. If the face contains boundary edges, these are removed. Situations
         may arise where the mesh is no longer manifold because the situation at a boundary vertex
         is not homeomorphic to a half disk. This, we can probably ignore since from the data
         structure point of view it is not really a problem that a vertex is incident on two holes - 
        a hole can be seen as a special type of face. The function returns false if the FaceID is 
         not valid, otherwise the function must complete. */
        bool remove_face(FaceID fid);

        /** Remove an edge from the Manifold.
            This function will remove the faces on either side and the edge itself in the process. Thus,
         it is a simple application of remove_face. */
        bool remove_edge(HalfEdgeID hid);

        /** Remove a vertex from the Manifold.
         This function merges all faces around the vertex into one and then removes 
         this resulting face. */
        bool remove_vertex(VertexID vid);
        
        /// number of  vertices
        size_t no_vertices() const { return kernel.no_vertices();}
        /// number of active faces
        size_t no_faces() const { return kernel.no_faces();}
        /// number of active halfedges
        size_t no_halfedges() const { return kernel.no_halfedges();}
        
        /// number of total vertices in kernel
        size_t allocated_vertices() const { return kernel.allocated_vertices();}
        /// number of total faces in kernel
        size_t allocated_faces() const { return kernel.allocated_faces();}
        /// number of total halfedges in kernel
        size_t allocated_halfedges() const { return kernel.allocated_halfedges();}
        
        /// check if ID of vertex is in use
        bool in_use(VertexID id) const { return kernel.in_use(id);}
        /// check if ID of face is in use
        bool in_use(FaceID id) const { return kernel.in_use(id);}
        /// check if ID of halfedge is in use
        bool in_use(HalfEdgeID id) const { return kernel.in_use(id);}
        
        IDIteratorPair<Vertex> vertices() const {return IDIteratorPair<Vertex>(kernel.vertices_begin(),
                                                                         kernel.vertices_end());}
        IDIteratorPair<Face> faces() const {return IDIteratorPair<Face>(kernel.faces_begin(),
                                                                         kernel.faces_end());}
        IDIteratorPair<HalfEdge> halfedges() const {return IDIteratorPair<HalfEdge>(kernel.halfedges_begin(),
                                                                         kernel.halfedges_end());}
        
        /// Iterator to first VertexID, optional argument defines if unused items should be skipped
        VertexIDIterator vertices_begin(bool skip = true) const { return kernel.vertices_begin();}
        /// Iterator to first FaceID, optional argument defines if unused items should be skipped
        FaceIDIterator faces_begin(bool skip = true) const { return kernel.faces_begin();}
        /// Iterator to first HalfEdgeID, optional argument defines if unused items should be skipped
        HalfEdgeIDIterator halfedges_begin(bool skip = true) const { return kernel.halfedges_begin();}
        
        /// Iterator to past the end VertexID
        VertexIDIterator vertices_end() const { return kernel.vertices_end();}
        /// Iterator topast the end FaceID
        FaceIDIterator faces_end() const { return kernel.faces_end();}
        /// Iterator to past the end HalfEdgeID
        HalfEdgeIDIterator halfedges_end() const {return kernel.halfedges_end(); }  

		
		/** \brief Bridge f0 and f1 by connecting the vertex pairs given in pairs.
		 This function creates a cylindrical connection between f0 and f1. f0 and f1 are removed and the vertices 
		 given in pairs are connected by edges. The result is a cylindrical connection that changes the genus of the object.
		 
		 This function leaves all error checking in the hands of the user (for now). The faces clearly should not have any
		 vertices or edges in common as this will create a non-manifold situation. Also the faces should face towards or away 
		 from each other and be in a position where it is reasonable to make the bridge. The connections should also make sense 
		 from a geometric point of view and should be in a counter clockwise loop on f0 and a clockwise loop on f1. No need to
		 connect all vertices.
		 
		 The function returns a vector of HalfEdgeIDs. Those are, of course, the connecting halfedges - also the opposite edges.
		 */
		std::vector<HalfEdgeID> bridge_faces(FaceID f0, FaceID f1, const std::vector<std::pair<VertexID, VertexID> >& pairs);


        /** \brief Collapse the halfedge h.
        The argument h is the halfedge being removed. The vertex v=h->opp->vert is the one being removed while h->vert survives.
        The final argument indicates whether the surviving vertex should have the average position of the former vertices.
        By default false meaning that the surviving vertex retains it position.
        This function is not guaranteed to keep the mesh sane unless, precond_collapse_edge has returned true !! */
        void collapse_edge(HalfEdgeID h, bool avg_vertices = false);

        /** \brief Split a face.
        The face, f, is split by creating an edge with endpoints v0 and v1 (the next two arguments). 
        The vertices of the old face between v0 and v1 (in counter clockwise order) continue to belong to f.
        The vertices between v1 and v0 belong to the new face. A handle to the new face is returned. */
        FaceID split_face_by_edge(FaceID f, VertexID v0, VertexID v1);

        /** \brief Split a polygon, f, by inserting a vertex at the barycenter.			
        This function is less likely to create flipped triangles than the split_face_triangulate function. 
        On the other hand, it introduces more vertices and probably makes the triangles more acute.
        A handle to the inserted vertex is returned. */
        VertexID split_face_by_vertex(FaceID f);
       // VertexID split_face_by_vertex(HalfEdgeID h);

        /** \brief Insert a new vertex on halfedge h.
        The new halfedge is insterted as the previous edge to h.
        A handle to the inserted vertex is returned. */
        VertexID split_edge(HalfEdgeID h);
        
        /** \brief Stitch two halfedges.
         Two boundary halfedges can be stitched together. This can be used to build a complex mesh
         from a bunch of simple faces. */
        bool stitch_boundary_edges(HalfEdgeID h0, HalfEdgeID h1);

        /** \brief Merges two faces into a single polygon. 
        The first face is f. The second face is adjacent to f along the halfedge h. 
        This function returns true if the merging was possible and false otherwise. 
        Currently merge only fails if the mesh is already illegal. Thus it should, in fact, never fail. */
        bool merge_faces(FaceID f, HalfEdgeID h);
		
		/** \brief Merge all faces in the one ring of a vertex into a single polygon.
		The vertex is given by v.
		 
		The return value is the FaceID of the resulting polygonal face. 
		InvalidFaceID is returned if 
		- the input vertex is not in use or 
		- the input vertex has valence less than two which is a degenerate case.
		- the input vertex is a boundary vertex of valence two - i.e. adjacent to just one face.
		- the same halfedge appears in two faces of the one ring of the input vertex: I.e.
		the input vertex is twice adjacent to the same face!
		 
		Note that this function can create some unusual and arguably degenerate meshes. For instance, 
		two triangles which share all vertices is collapsed to a single pair of vertices connected by 
		a pair of halfedges bounding the same face. */
		FaceID merge_one_ring(VertexID v, float max_loop_length = FLT_MAX);

        /** \brief Close hole given by the invalid face of halfedgehandle h.
         returns FaceID of the created face or the face that is already there if the 
         face was not InvalidFaceID. */
        FaceID close_hole(HalfEdgeID h);
        
        /** \brief Create hole by opening a slit at a vertex.
        This function creates a hole by slitting the mesh at a vertex v along a pair of halfedges.
        The function returns the vertex created by this operation.
        The first halfedge h_in is oriented towards the vertex and the second, h_out, away from the 
        vertex. Neither h_in nor h_out may be boundary halfedges for a hole since the result would be
        two holes separated by an edge. That is an abomination we should avoid. In particular because
        it would be invisible. h_in and h_out also should not be each other's opposite edges. That would 
         result in an isolated vertex. */
         VertexID slit_vertex(VertexID v, HalfEdgeID h_in, HalfEdgeID h_out);
        
        /** \brief Create hole by slitting open the mesh along the path given as argument.
         This function returns the HalfEdgeID of one of the halfedges bounding the created hole.
         The path is specified as a selection set of vertices. If the selected vertices form a closed 
         loop, a piece is cut off from the mesh. Surprising results can occur if the selected vertices
         can be connected by more than one sequence (or a self intersecting sequence) of edges. */
        HalfEdgeID slit_edges(VertexAttributeVector<int>& selection);

        /// \brief Flip an edge h. 
        void flip_edge(HalfEdgeID h);

        /// Return reference to position given by VertexID
        Vec& pos(VertexID id);
        /// Return const reference to position given by VertexID
        const Vec& pos(VertexID id) const;
        
        /// Return a reference to the entire positions attribute vector
        VertexAttributeVector<Vec>& positions_attribute_vector();

        /// Return a const reference to the entire positions attribute vector
        const VertexAttributeVector<Vec>& positions_attribute_vector() const;
        
        /// Clear the mesh. Remove all faces, halfedges, and vertices.
        void clear();

        /// Remove unused items from Mesh, map argument is to be used for attribute vector cleanups in order to maintain sync.
        void cleanup(IDRemap& map);
        /// Remove unused items from Mesh
        void cleanup();
        
        /// Returns a Walker to the out halfedge of vertex given by VertexID
        Walker walker(VertexID id) const;
        /// Returns a Walker to the last halfedge of face given by FaceID
        Walker walker(FaceID id) const;
        /// Returns a Walker to the halfedge given by HalfEdgeID
        Walker walker(HalfEdgeID id) const;
        

    private:
        
        ConnectivityKernel kernel;
        
        VertexAttributeVector<Vec> positions;

        // private template for building the manifold from various types
        template<typename size_type, typename float_type, typename int_type>
        void build_template(size_type no_vertices,
                            const float_type* vertvec,
                            size_type no_faces,
                            const int_type* facevec,
                            const int_type* indices);

        /// Set the next and prev indices of the first and second argument respectively.
        void link(HalfEdgeID h0, HalfEdgeID h1);

        /// Glue halfedges by letting the opp indices point to each other.
        void glue(HalfEdgeID h0, HalfEdgeID h1);

        /// Auxiliary function called from collapse
        void remove_face_if_degenerate(HalfEdgeID h);

        /// Ensure boundary consistency.
        void ensure_boundary_consistency(VertexID v);
    };

    /** \brief Verify Manifold Integrity
    Performs a series of tests to check that this is a valid manifold.
    This function is not rigorously constructed but seems to catch all problems so far. 
    The function returns true if the mesh is valid and false otherwise. */
    bool valid(const Manifold& m);

    /// Calculate the bounding box of the manifold
    void bbox(const Manifold& m, Manifold::Vec& pmin, Manifold::Vec& pmax);

    /// Calculate the bounding sphere of the manifold
    void bsphere(const Manifold& m, Manifold::Vec& c, float& r);

    /** \brief Test for legal edge collapse.
    The argument h is the halfedge we want to collapse. 
    If this function does not return true, it is illegal to collapse h. 
    The reason is that the collapse would violate the manifold property of the mesh.
    The test is as follows:
    1.  For the two vertices adjacent to the edge, we generate a list of all their neighbouring vertices. 
    We then generate a  list of the vertices that occur in both these lists. 
    That is, we find all vertices connected by edges to both endpoints of the edge and store these in a list.
    2.  For both faces incident on the edge, check whether they are triangular.
    If this is the case, the face will be removed, and it is ok that the the third vertex is connected to both endpoints. 
    Thus the third vertex in such a face is removed from the list generated in 1.
    3.  If the list is now empty, all is well. 
    Otherwise, there would be a vertex in the new mesh with two edges connecting it to the same vertex. Return false.
    4.  TETRAHEDRON TEST:
    If the valency of both vertices is three, and the incident faces are triangles, we also disallow the operation. 
    Reason: A vertex valency of two and two triangles incident on the adjacent vertices makes the construction collapse.
    5.  VALENCY 4 TEST:
    If a triangle is adjacent to the edge being collapsed, it disappears.
    This means the valency of the remaining edge vertex is decreased by one.
    A valency two vertex reduced to a valency one vertex is considered illegal.
    6.  PREVENT MERGING HOLES:
    Collapsing an edge with boundary endpoints and valid faces results in the creation where two holes meet.
    A non manifold situation. We could relax this...
	7. New test: if the same face is in the one-ring of both vertices but not adjacent to the common edge,
	then the result of a collapse would be a one ring where the same face occurs twice. This is disallowed as the resulting
	 face would be non-simple.	*/
    bool precond_collapse_edge(const Manifold& m, HalfEdgeID h);

    /** \brief Test fpr legal edge flip. 
    Returns false if flipping cannot be performed. This is due to one of following: 
    1.  one of the two adjacent faces is not a triangle. 
    2.  Either end point has valency three.
    3.  The vertices that will be connected already are. */
    bool precond_flip_edge(const Manifold& m, HalfEdgeID h);

    /// Returns true if the halfedge is a boundary halfedge.
    bool boundary(const Manifold& m, HalfEdgeID h);

    /// Return the geometric length of a halfedge.
    double length(const Manifold& m, HalfEdgeID h);

    /// Returns true if the vertex is a boundary vertex.
    bool boundary(const Manifold& m, VertexID v);

    /// Compute valency, i.e. number of incident edges.
    int valency(const Manifold& m, VertexID v);

    /// Compute the vertex normal. This function computes the angle weighted sum of incident face normals.
    Manifold::Vec normal(const Manifold& m, VertexID v);

    /// Returns true if the two argument vertices are in each other's one-rings.
    bool connected(const Manifold& m, VertexID v0, VertexID v1);

    /// Compute the number of edges of a face
    int no_edges(const Manifold& m, FaceID f);

    /** Compute the normal of a face. If the face is not a triangle,
    the normal is not defined, but computed using the first three
    vertices of the face. */
    Manifold::Vec normal(const Manifold& m, FaceID f);

    /// Compute the area of a face. 
    double area(const Manifold& m, FaceID f);

    /// Compute the perimeter of a face. 
    double perimeter(const Manifold& m, FaceID f);

    /// Compute the centre of a face
    Manifold::Vec centre(const Manifold& m, FaceID f);

    /*******************************************************************
    * Manifold code
    *******************************************************************/

    inline Manifold::Manifold(){}

    inline Manifold::Vec& Manifold::pos(VertexID id)
    { return positions[id]; }
    inline const Manifold::Vec& Manifold::pos(VertexID id) const
    { return positions[id]; }
    
    inline VertexAttributeVector<Manifold::Vec>& Manifold::positions_attribute_vector()
    {
        return positions;
    }
    
    inline const VertexAttributeVector<Manifold::Vec>& Manifold::positions_attribute_vector() const
    {
        return positions;    
    }

    inline void Manifold::clear()
    { 
        kernel.clear();
        positions.clear();
    }

    inline Walker Manifold::walker(VertexID id) const
    { return Walker(kernel, kernel.out(id)); }
    inline Walker Manifold::walker(FaceID id) const
    { return Walker(kernel, kernel.last(id)); }
    inline Walker Manifold::walker(HalfEdgeID id) const
    { return Walker(kernel, id); }


    inline void Manifold::cleanup(IDRemap& map)
    {   
        kernel.cleanup(map);
        positions.cleanup(map.vmap);
    }
    
    inline void Manifold::cleanup()
    {
        IDRemap map;
        Manifold::cleanup(map);
    }
    
    inline int circulate_vertex_ccw(const Manifold& m, VertexID v, std::function<void(Walker&)> f)
    {
        Walker w = m.walker(v);
        for(; !w.full_circle(); w = w.circulate_vertex_ccw()) f(w);
        return w.no_steps();
    }
    inline int circulate_vertex_ccw(const Manifold& m, VertexID v, std::function<void(VertexID)> f)
    {
        return circulate_vertex_ccw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.vertex());});
    }
    inline int circulate_vertex_ccw(const Manifold& m, VertexID v, std::function<void(FaceID)> f)
    {
        return circulate_vertex_ccw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.face());});
    }
    inline int circulate_vertex_ccw(const Manifold& m, VertexID v, std::function<void(HalfEdgeID)> f)
    {
        return circulate_vertex_ccw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.halfedge());});
    }
    
    inline int circulate_vertex_cw(const Manifold& m, VertexID v, std::function<void(Walker&)> f)
    {
        Walker w = m.walker(v);
        for(; !w.full_circle(); w = w.circulate_vertex_cw()) f(w);
        return w.no_steps();
    }
    inline int circulate_vertex_cw(const Manifold& m, VertexID v, std::function<void(VertexID)> f)
    {
        return circulate_vertex_cw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.vertex());});
    }
    inline int circulate_vertex_cw(const Manifold& m, VertexID v, std::function<void(FaceID)> f)
    {
        return circulate_vertex_cw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.face());});
    }
    inline int circulate_vertex_cw(const Manifold& m, VertexID v, std::function<void(HalfEdgeID)> f)
    {
        return circulate_vertex_cw(m, v, (std::function<void(Walker&)>)[&](Walker& w){f(w.halfedge());});
    }
    
    inline int circulate_face_ccw(const Manifold& m, FaceID f, std::function<void(Walker&)> g)
    {
        Walker w = m.walker(f);
        for(; !w.full_circle(); w = w.circulate_face_ccw()) g(w);
        return w.no_steps();
    }
    inline int circulate_face_ccw(const Manifold& m, FaceID f, std::function<void(VertexID)> g)
    {
        return circulate_face_ccw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.vertex());});
    }
    inline int circulate_face_ccw(const Manifold& m, FaceID f, std::function<void(FaceID)> g)
    {
        return circulate_face_ccw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.opp().face());});
    }
    inline int circulate_face_ccw(const Manifold& m, FaceID f, std::function<void(HalfEdgeID)> g)
    {
        return circulate_face_ccw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.halfedge());});
    }
    
    inline int circulate_face_cw(const Manifold& m, FaceID f, std::function<void(Walker&)> g)
    {
        Walker w = m.walker(f);
        for(; !w.full_circle(); w = w.circulate_face_cw()) g(w);
        return w.no_steps();
    }
    inline int circulate_face_cw(const Manifold& m, FaceID f, std::function<void(VertexID)> g)
    {
        return circulate_face_cw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.vertex());});
    }
    inline int circulate_face_cw(const Manifold& m, FaceID f, std::function<void(FaceID)> g)
    {
        return circulate_face_cw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.opp().face());});
    }
    inline int circulate_face_cw(const Manifold& m, FaceID f, std::function<void(HalfEdgeID)> g)
    {
        return circulate_face_cw(m, f, (std::function<void(Walker&)>)[&](Walker& w){g(w.halfedge());});
    }
    

}
