/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */
/**
 * @file ManifoldRenderer.h
 * @brief Tool for rendering a HMesh Manifold.
 */

#ifndef __MESHEDIT_RENDERER_H__
#define __MESHEDIT_RENDERER_H__

#include "../CGLA/Vec2d.h"
#include "../GL/glew.h"
#include "../HMesh/harmonics.h"
#include "../GLGraphics/draw.h"
#include "../GLGraphics/Console.h"
#include "../GLGraphics/IDBufferWireFrameRenderer.h"
#include "../CGLA/Vec4d.h"

namespace HMesh
{
    template<typename ITEM>
    class VertexAttributeVector;
}

namespace GLGraphics {
    
    /** Ancestral class for Manifold rendering. Do not use directly. Its only purpose is to
     create a display list and remove it when the object is destroyed. This is an example
     of the RAII "resource acquisition is initialization" idiom */
    class ManifoldRenderer
	{
	protected:
		GLuint display_list;
	public:
		ManifoldRenderer(): display_list(glGenLists(1))	{}
		virtual ~ManifoldRenderer()
		{
			glDeleteLists(display_list, 1);
		}
        /// Produce a display list containing geometry and normals (which may be smooth or per face).
		virtual void compile_display_list(const HMesh::Manifold& m, bool smooth) {}

		virtual void draw()
		{
			glCallList(display_list);
		}
	};
    
    
    /** Wireframe rendering. This is a nasty complex class that relies on other classes. The trouble
     is that for non-triangle meshes, we need to use another approach than for triangle meshes.
     This class is really a front end for a couple of other classes. */
    class WireframeRenderer: public ManifoldRenderer
	{
		GLGraphics::IDBufferWireframeRenderer* idbuff_renderer;
		
		int maximum_face_valency(const HMesh::Manifold& m);
		
	public:
		~WireframeRenderer()
		{
			delete idbuff_renderer;
		}
		
		WireframeRenderer(HMesh::Manifold& m, bool flat);
		
		void draw();
	};
    
    /** SimpleShaderRenderer is a very basic class for drawing a Manifold with shading.
     It is a convenient way to draw a surface using vertex and fragment shaders since it takes
     care of initializing the shaders and producing a display list for the geometry.
     
     Geometry shaders typically add more complexity and are left out of this class, so you cannot add a
     geometry shader.
     
     While this class can be used directly, the normal procedure is to inherit from SimpleShaderRenderer
     and then pass the shaders to the constructor. The strings defining the shaders would fit
     nicely as static constant strings (see e.g. ToonRenderer or GlazedRenderer) in your inherited class.
     
     If you need to define more attributes or uniforms, you need to take charge. Your inherited class's
     constructor should then use the default constructor of SimpleShaderRenderer. You can call init_shaders
     to initialize the shaders and then compile the display list yourself with the needed uniforms and
     attributes - rather than calling compile_display_list which only puts vertices and normals in the list.
     */
    class SimpleShaderRenderer: public ManifoldRenderer
	{
		/// Compile the vertex and fragment shaders and link to form shader program.
		void init_shaders(const std::string& vss,
						  const std::string& fss);
		
        
	protected:
        
		GLuint prog,vs,fs;
        
	public:
				
		/** This constructor simply initializes the shaders. It does not create the display list.
         Use if you shader has extra attributes. */
		SimpleShaderRenderer(const std::string& vss,
							 const std::string& fss) {init_shaders(vss,fss);}
		
		/// Produce a display list containing geometry and normals (which may be smooth or per face).
		virtual void compile_display_list(const HMesh::Manifold& m, bool smooth);

		/// Releases the program and shaders.
		~SimpleShaderRenderer()
		{
			glDeleteProgram(prog);
			glDeleteShader(vs);
			glDeleteShader(fs);
		}
		
		/// Do the actual drawing. Simply calls the display list if this function is not overloaded.
		virtual void draw();
        
	};
    
    /** Ugly basic gouraud rendering. This class uses OpenGL's fixed function pipeline. */
    class NormalRenderer: public SimpleShaderRenderer
    {
        const static std::string vss;
        const static std::string fss;
        
    public:
        NormalRenderer():
        SimpleShaderRenderer(vss, fss) {}
    };
    
    /** Debug renderer. Color code faces and show vertices as balls.*/
    class DebugRenderer: public SimpleShaderRenderer
    {
        const static std::string vss;
        const static std::string fss;
        
    public:
        static HMesh::VertexAttributeVector<CGLA::Vec3f> vertex_colors;
        static HMesh::FaceAttributeVector<CGLA::Vec3f> face_colors;
        static HMesh::HalfEdgeAttributeVector<CGLA::Vec3f> edge_colors;
    public:
        DebugRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m, bool smooth, float rad);
        
    };
    
    /** Render reflection lines. This class renders the object as if it is specular and inside
     an infinitely long, vertical cylinder with white strips (also vertical). Useful if you
     want to see whether the surface is smooth or has kinks. */
    class ReflectionLineRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		ReflectionLineRenderer(): SimpleShaderRenderer(vss, fss) {}

	};
    
    /** Render isophotes with respect to a lightsource in the eye. Useful if you
     want to see whether the surface is smooth or has kinks. */
    class IsophoteLineRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		IsophoteLineRenderer(): SimpleShaderRenderer(vss, fss) {}
		
	};
    
    /** The toon renderer simply quantizes the shading to give a toonish appearance
     with a fat black silhouette. */
    
    class ToonRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		ToonRenderer(): SimpleShaderRenderer(vss, fss) {}
		
	};
    
    /** Render like glazed ceramics. Looks cool. I will add more to this. */
    class GlazedRenderer: public SimpleShaderRenderer
	{
		float bsphere_rad;
		const static std::string vss;
		const static std::string fss;
	public:
		GlazedRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m, bool smooth);
	};
    
    /** Render a scalar field. Positive scalars are mapped to blue and negative to red.
     This class also has controls for gamma correction which is highly useful if the
     scalars are mostly small or large and simply scaling to the 0-1 range does not
     produce a good result. */
    class ScalarFieldRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		ScalarFieldRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m, bool smooth,
                                  HMesh::VertexAttributeVector<double>& field,
                                  double min_val,
                                  double max_val,
                                  float gamma = 2.2,
                                  int use_stripes=0,
                                  int color_signed=0,
                                  int use_shading=0);
	};
    
    /** Render a checkerboard pattern based on input texture map */
    class CheckerBoardRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
        static HMesh::VertexAttributeVector<CGLA::Vec2f> param;
		CheckerBoardRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m, bool smooth);
	};

    
    
    /** Ambient occlusion renderer. Very similar to ScalarFieldRender. Simply assumes that the input values are
     mean curvatures which in some sense indicate how concave the surface is.*/
    class AmbientOcclusionRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		AmbientOcclusionRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& field, double max_val);
	};
    
    
    /** Line fields are rendered by convolving a noise function in the direction of the line.
     This is useful, for instance, for curvature rendering. */
    class LineFieldRenderer: public SimpleShaderRenderer
	{
		const static std::string vss;
		const static std::string fss;
	public:
		LineFieldRenderer(): SimpleShaderRenderer(vss, fss) {}
        void compile_display_list(const HMesh::Manifold& m,HMesh::VertexAttributeVector<CGLA::Vec3d>& lines);
	};
    
    class HarmonicsRenderer: public GLGraphics::ManifoldRenderer
    {
        static GLuint prog_P0;
        static GLGraphics::Console::variable<float> display_harmonics_time;
        static GLGraphics::Console::variable<int> display_harmonics_diffuse;
        static GLGraphics::Console::variable<int> display_harmonics_highlight;
        static GLGraphics::Console::variable<int> display_harmonics_e0;
        static GLGraphics::Console::variable<int> display_harmonics_e1;
        
        HMesh::Manifold* m;
        HMesh::Harmonics* h;
        
        /// Draw with eigenvalues
        void draw_adf();
        void draw_esum();
        
    public:
        HarmonicsRenderer(HMesh::Manifold& _m, HMesh::Harmonics* _h, GLGraphics::Console& cs);
        
        /// Parse keystrokes that would influence the interactive display
        void parse_key(unsigned char key);
        
    };

    
}

#endif