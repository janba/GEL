/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file IDBufferWireFrameRenderer.h
 * @brief Render a Manifold in wireframe. Not fast - but general.
 */

#ifndef __GLGRAPHICS_IDBUFFERWIREFRAMERENDERER_H__
#define __GLGRAPHICS_IDBUFFERWIREFRAMERENDERER_H__

#include "../GL/glew.h"
#include "../CGLA/Vec3uc.h"
#include "../CGLA/Vec3f.h"
#include "../HMesh/Manifold.h"

namespace GLGraphics
{
	/// Class for ID buffer based wireframe rendering. Handles all types of gons. Use with Manifold
	class IDBufferWireframeRenderer
		{
			HMesh::Manifold* mesh;
			
			GLint id_attrib;
			GLint popp_attrib;
			GLint disp_attrib;
			int XSZ, YSZ;
			float thickness;
			float transition;
			
			GLuint vs,fs;
			GLuint line_prog;
			GLuint idmap;
			
			GLuint vertex_buffername;
			GLuint colors_buffername;			
			GLuint line_id_attrib;
			GLuint line_vertex_pos;
			GLuint line_disp_attrib;
			GLuint line_opp_attrib;
			
			int triangles, quads;
		public:
			
			CGLA::Vec3uc id_get(unsigned int i)
			{
				i = i+1;
				return CGLA::Vec3uc(i&0xff, (i&0xff00)/256, (i&0xff0000)/65536);
			}
			
			IDBufferWireframeRenderer(int _XSZ, int _YSZ, 
							  HMesh::Manifold& mesh,
							  float _thickness=0.0, 
							  float _transition=1.8,
							  int atten_mode=0);
			
			~IDBufferWireframeRenderer();
			
			void draw(const CGLA::Vec3f& color, 
					const CGLA::Vec3f& clear_color);
		};
	
}
#endif
