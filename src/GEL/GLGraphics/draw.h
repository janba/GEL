/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file draw.h
 * @brief Draw various GEL entities.
 */

#ifndef __GLGRAPHICS_DRAW_H_
#define __GLGRAPHICS_DRAW_H_

#include "../Geometry/TriMesh.h"
#include "../Geometry/AABox.h"
#include "../Geometry/OBox.h"
#include "../Geometry/BoundingINode.h"
#include "../Geometry/BoundingLNode.h"
#include "../Geometry/BoundingTree.h"

namespace HMesh
{
    class Manifold;
}

namespace GLGraphics
{
    /// Simply draw a tessellated ball.
    void draw_ball();
    
	/// Draw an indexed face set (just a triangle list) with normals computed on the fly.
	void draw(const Geometry::IndexedFaceSet& geometry);

	/// Draw a triangles mesh. Inefficient function that should be compiled into a display list.
	void draw(const Geometry::TriMesh& tm, bool per_vertex_norms=true);
	
	/// Load textures if available
	void load_textures(Geometry::TriMesh& tm);	
	
	/// Draw a HMesh Manifold. Inefficient and should be compiled into display list
	void draw(const HMesh::Manifold& m, bool per_vertex_norms=true);

	
	/// Draw an axis aligned bounding box
	void draw(const Geometry::AABox& box);
	
	/// Draw an oriented bounding box
	void draw(const Geometry::OBox& box);
	
	/// Draw an object of type T which contains only triangles as wireframe. In practice T = Manifold or TriMesh.
	template<typename T>
	inline void draw_triangles_in_wireframe(T& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);

  /// Draw an object of type T as wireframe.  In practice T = Manifold or TriMesh.
	template<class T>
	void draw_wireframe_oldfashioned(const T& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);
	
	
	template<class BoxType>
	void draw(const Geometry::BoundingINode<BoxType>& node, int level, int max_level);
	template<class BoxType>
    void draw(const Geometry::BoundingLNode<BoxType>& node, int level, int max_level);
	template<class BoxType>
    void draw(const Geometry::BoundingTree<BoxType>& tree, int max_level = 1e6);

    /// Read depth buffer and extract the depth of a pixel
    bool depth_pick(int x, int y, float& depth);
    
    /** Convert depth value to 3D point. Supposes that the modelview and projection matrices
     are set in the same way as when rendering. */
    CGLA::Vec3d screen2world(int x, int y, float depth);
    
    /** Project the point given as argument using modelview and projection matrices */
    CGLA::Vec3d world2screen(const CGLA::Vec3d& p);


}
#endif
