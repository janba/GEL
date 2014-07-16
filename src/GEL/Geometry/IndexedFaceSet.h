/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file IndexedFaceSet.h
 * @brief Simple indexed triangle mesh data structure.
 */

#ifndef __GEOMETRY_GEOMETRY_INDEXED_FACE_SET_H__
#define __GEOMETRY_GEOMETRY_INDEXED_FACE_SET_H__

#include <vector>
#include "../CGLA/Vec3f.h"
#include "../CGLA/Vec3i.h"

namespace Geometry
{
	const CGLA::Vec3i NULL_FACE(-1,-1,-1);
	
	/** \brief This class represents the simplest possible triangle mesh data
			structure.

			Faces must be triangular, and 3D coordinates are the only values 
			associated with the vertices. */
	class IndexedFaceSet
	{
		/// Vertices
    std::vector<CGLA::Vec3f> verts;
			
		/// Vec3is (which must be triangles)
    std::vector<CGLA::Vec3i> faces;
				
	public:

		IndexedFaceSet(): verts(0), faces(0) {}

		// ----------------------------------------
		// Functions that operate on faces
		// ----------------------------------------
				
		/** Add a face and return the index of the new face. 
				
		If the optional idx argument is present, the face array is resized so 
		that the new index == idx. */
		int add_face(const CGLA::Vec3i& f, int idx=-1) 
		{
			if(idx < 0)
				idx = static_cast<int>(faces.size());
			faces.resize(idx+1,NULL_FACE);
			faces[idx] = f;
			return idx;
		}

		/// Return the number of faces.
		int no_faces() const {return static_cast<int>(faces.size());}

		/** Return the face corresponding to a given index. 
				
				The NULL_FACE	is returned if the index is out of bounds. */
		const CGLA::Vec3i& face(size_t idx) const
		{
			if(idx<faces.size())
				return faces[idx];
			return NULL_FACE;
		}

		/// Assign f to a face of index idx
		CGLA::Vec3i& face_rw(int idx) 
		{
			return faces[idx];
		}

		// ----------------------------------------
		// Functions that operate on vertices
		// ----------------------------------------

		/// Add a vertex and return the index of the vertex.
		int add_vertex(const CGLA::Vec3f& v)
		{
			int idx= static_cast<int>(verts.size());
			verts.push_back(v);
			return idx;
		}

		/// Return the number of vertices.
		int no_vertices() const {return static_cast<int>(verts.size());}

		/** Return the vertex corresponding to a given index. 
				
		User is responsible	for bounds checking. */
		const CGLA::Vec3f& vertex(int idx) const
		{
			return verts[idx];
		}

		/// Assign v to a vertex of index idx
		CGLA::Vec3f& vertex_rw(int idx)
		{
			return verts[idx];
		}

	};
}

#endif
