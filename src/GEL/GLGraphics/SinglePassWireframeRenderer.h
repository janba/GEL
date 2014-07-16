/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file SinglePassWireframeRenderer.h
 * @brief Render triangles in wire frame using efficient single pass method.
 */

#ifndef __GLGRAPHICS_SINGLE_PASS_WIREFRAME__
#define __GLGRAPHICS_SINGLE_PASS_WIREFRAME__

#include "../GL/glew.h"
#include "../CGLA/Vec3f.h"

namespace GLGraphics
{
	/** Class for wireframe rendering. Enable it will case all triangles to be drawn to be drawn as
		wireframe. Only triangles are supported, but it is fast. */
	class SinglePassWireframeRenderer
		{
			GLhandleARB prog, old_prog;
		public:
			/// The constructor creates the wireframe shader program. It is static so the program is shared.
			SinglePassWireframeRenderer();
			
			/// Enable the wireframe shader with an optional line color
			bool enable(const CGLA::Vec3f& line_color = CGLA::Vec3f(1,0,0));
			
			/// Switch back to the shader program which was enabled at the time of calling enable.
			void disable();
		};
}
#endif