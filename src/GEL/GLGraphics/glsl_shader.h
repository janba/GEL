/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file glsl_shader.h
 * @brief Simple functions for loading GLSL shaders.
 */

#ifndef __GLGRAPHICS_GLSL_SHADER_H__
#define __GLGRAPHICS_GLSL_SHADER_H__

#include "../GL/glew.h"
#include <string>

/* It is a little tricky to make shader programs in C++ using GLSL but the problems are
  almost all silly little things. For instance, how do you robustly read from a text 
  file? How do you compile a shader and check for errors? In OpenGL What is the difference 
  between a GLSL program and a shader anyway?

  The short answer to the last question is this: A "shader" can be either a vertex shader, 
  a geometry shader (in OpenGL 2.0) or a fragment shader. A "program" is a linked combination
  of these three types of shaders. Note that you don't have to have a geometry shader.
  
  This API attempts to obviate the need for an answer to the first two questions by providing 
  an API for loading and compiling shaders and checking for errors. However, I don't
  include functions for creating the program, attaching shaders and linking. Why not??!
  
  Because the need for flexibility means that the API would be just as complex as just using
  the OpenGL functions directly! However, loading shaders and checking for errors in compiled shaders
  is different. It makes sense to wrap that. It also makes sense to wrap the error checking for 
  programs that are linked, so there is a function for that too.
  
  Since shader loading and error check sandwhich the two calls needed for compilation, the most 
  important function in this API, create_glsl_shader, loads a shader, compiles it, checks for errors 
  and returns the shader handle. There is also a version which creates a shader from a string.
  
  There is some code below to illustrate usage.
  
<code snippet>
  // Create shaders directly from file
  GLuint vs = create_glsl_shader(GL_VERTEX_SHADER, shader_path, "tri.vert");
  GLuint gs = create_glsl_shader(GL_GEOMETRY_SHADER_EXT, shader_path, "tri.geom");
  GLuint fs = create_glsl_shader(GL_FRAGMENT_SHADER, shader_path, "tri.frag");

  // Create the program
  prog_P0 = glCreateProgram();
  
  // Attach all shaders
  if(vs) glAttachShader(prog_P0, vs);
  if(gs) glAttachShader(prog_P0, gs);
  if(fs) glAttachShader(prog_P0, fs);
  
  // Specify input and output for the geometry shader. Note that this must be
  // done before linking the program.
  glProgramParameteriEXT(prog_P0,GL_GEOMETRY_INPUT_TYPE_EXT,GL_TRIANGLES);
  glProgramParameteriEXT(prog_P0,GL_GEOMETRY_VERTICES_OUT_EXT,3);
  glProgramParameteriEXT(prog_P0,GL_GEOMETRY_OUTPUT_TYPE_EXT,GL_TRIANGLE_STRIP);

  // Link the program object and print out the info log
  glLinkProgram(prog_P0);
  print_glsl_program_log(prog_P0);
  
  // Install program object as part of current state
  glUseProgram(prog_P0);

  // Set the value of a uniform
  glUniform2f(glGetUniformLocation(prog_P0,"WIN_SCALE"), win_size_x/2.0, win_size_y/2.0);
</code snippet>
  
  Happy shader coding.
  
  Andreas B¾rentzen, 2007
  
  */

namespace GLGraphics
{
  /// Print the info log for a program if the status is not OK.
  void print_glsl_program_log(GLuint program);

  /// Print the info log for a shader if the status is not OK.
  void print_glsl_shader_log(GLuint shader);
  
  /** The two arguments are concatenated to form the name with full path of a text file.
    This file is read and returned as a string. */
  const std::string read_glsl_source(const std::string& path, const std::string& file);
  
  /** Create a shader of type specified by the first argument from a source string given
    as second argument.  Return shader handle. If there is a problem, the infolog is 
    printed and 0 is returned (unless the third argument is false). */
  GLuint create_glsl_shader(GLuint stype, const std::string& src, bool print_log = true);
    
  /** Create a shader of type specified as first argument from the file indicated by the 
    supplied path and file name (second and third arguments) and return a shader handle. 
    This function is only a convenience function wrapping read_glsl_source and 
    create_glsl_shader */
  GLuint create_glsl_shader(GLuint stype, const std::string& path, const std::string& file);
}

#endif
