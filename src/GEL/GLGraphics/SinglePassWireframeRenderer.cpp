/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>
#include "../GLGraphics/glsl_shader.h"
#include "SinglePassWireframeRenderer.h"

using namespace CGLA;
using namespace std;
using namespace GLGraphics;

namespace 
{
	const string vp = 
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"varying vec4 diffuseIn;\n"
	"void main(void)\n"
	"{\n"
	"   float ndot = dot(gl_NormalMatrix*gl_Normal,vec3(0,0,1));\n"
	"   diffuseIn = vec4(0.7,0.9,1,1) * abs(ndot);\n"
	"   gl_Position =  ftransform();\n"
	"}\n";
	
	const string gp = 
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"#extension GL_EXT_geometry_shader4 : enable\n"
	"\n"
	"uniform vec2 WIN_SCALE;\n"
	"varying in vec4 diffuseIn[3];\n"
	"varying vec4 diffuse;\n"
	"noperspective varying vec3 dist;\n"
	"void main(void)\n"
	"{\n"
	"  vec2 p0 = WIN_SCALE * gl_PositionIn[0].xy/gl_PositionIn[0].w;\n"
	"  vec2 p1 = WIN_SCALE * gl_PositionIn[1].xy/gl_PositionIn[1].w;\n"
	"  vec2 p2 = WIN_SCALE * gl_PositionIn[2].xy/gl_PositionIn[2].w;\n"
	"  \n"
	"  vec2 v0 = p2-p1;\n"
	"  vec2 v1 = p2-p0;\n"
	"  vec2 v2 = p1-p0;\n"
	"  float area = abs(v1.x*v2.y - v1.y * v2.x);\n"
	"\n"
	"  dist = vec3(area/length(v0),0,0);\n"
	"  gl_Position = gl_PositionIn[0];\n"
	"  diffuse = diffuseIn[0];\n"
	"  EmitVertex();\n"
	"	\n"
	"  dist = vec3(0,area/length(v1),0);\n"
	"  gl_Position = gl_PositionIn[1];\n"
	"  diffuse = diffuseIn[1];\n"
	"  EmitVertex();\n"
	"\n"
	"  dist = vec3(0,0,area/length(v2));\n"
	"  gl_Position = gl_PositionIn[2];\n"
	"  diffuse = diffuseIn[2];\n"
	"  EmitVertex();\n"
	"\n"
	"  EndPrimitive();\n"
	"}\n";
	
	const string fp =
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"\n"
	"noperspective varying vec3 dist;\n"
	"varying vec4 diffuse;\n"
	"uniform vec4 WIRE_COL;\n"
	"\n"
	"void main(void)\n"
	"{\n"
	"	float d = min(dist[0],min(dist[1],dist[2]));\n"
	"	// Compute intensity\n"
	"\n"
	" 	float I = exp2(-2*d*d);\n"
	" 	gl_FragColor = I*WIRE_COL + (1.0 - I)*diffuse;\n"
	"}\n";
	
}


namespace GLGraphics
{
	
	SinglePassWireframeRenderer::SinglePassWireframeRenderer()
	{		
		static bool was_here = false;
		if(!was_here)
		{
			was_here = true;
			// Create s	haders directly from file
			static GLuint vs = create_glsl_shader(GL_VERTEX_SHADER, vp);
			static GLuint gs = create_glsl_shader(GL_GEOMETRY_SHADER_EXT, gp);
			static GLuint fs = create_glsl_shader(GL_FRAGMENT_SHADER, fp);
			
			// Create the program
			static GLuint _prog = glCreateProgram();
			prog = _prog;
			
			// Attach all shaders
			if(vs) glAttachShader(prog, vs);
			if(gs) glAttachShader(prog, gs);
			if(fs) glAttachShader(prog, fs);
			
			// Specify input and output for the geometry shader. Note that this must be
			// done before linking the program.
			glProgramParameteriEXT(prog,GL_GEOMETRY_INPUT_TYPE_EXT,GL_TRIANGLES);
			glProgramParameteriEXT(prog,GL_GEOMETRY_VERTICES_OUT_EXT,3);
			glProgramParameteriEXT(prog,GL_GEOMETRY_OUTPUT_TYPE_EXT,GL_TRIANGLE_STRIP);
			
			// Link the program object and print out the info log
			glLinkProgram(prog);
		}
	}
	
	bool SinglePassWireframeRenderer::enable(const CGLA::Vec3f& line_col)
	{
		GLint o;
		glGetIntegerv(GL_CURRENT_PROGRAM, &o);
		old_prog = o;
		glUseProgram(prog);
		
		GLint vpdim[4];
		glGetIntegerv(GL_VIEWPORT,vpdim); 
		
		// Set the value of a uniform
		glUniform2f(glGetUniformLocation(prog,"WIN_SCALE"), 
					static_cast<float>(vpdim[2]/2), static_cast<float>(vpdim[3]/2));
		glUniform4f(glGetUniformLocation(prog,"WIRE_COL"), line_col[0], line_col[1], line_col[2], 0);
		
		return true;
	}
	
	void SinglePassWireframeRenderer::disable()
	{
		glUseProgram(old_prog);
	}
	
}
