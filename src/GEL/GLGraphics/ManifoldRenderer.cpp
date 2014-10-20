/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ManifoldRenderer.h"

#if !defined (WIN32)
  #include <GLUT/GLUT.h>
#else
  #include <GL/GLUT.h>
#endif


#include <algorithm>
#include <string>
#include <cstdlib>
#include "../Geometry/TriMesh.h"
#include "../CGLA/Mat3x3d.h"
#include "../GLGraphics/glsl_shader.h"
#include "../GLGraphics/draw.h"
#include "../HMesh/Manifold.h"
#include "../HMesh/AttributeVector.h"
#include "../HMesh/curvature.h"

using namespace CGLA;
using namespace HMesh;
using namespace std;
using namespace Geometry;
namespace GLGraphics
{    
    GLuint get_noise_texture_id()
    {
        static GLuint texname=0;
        static bool was_here = false;
        
        if(!was_here)
        {
            was_here = true;
            int width = 32;
            int height = 32;
            int depth = 32;
            vector<unsigned char> texels(width*height*depth);
            for (int i = 0; i < width*height*depth; ++i)
            {
                int intensity = 255.0 * (float(gel_rand()) / GEL_RAND_MAX);
                texels[i] = (unsigned char) intensity;
            }
            
            glGenTextures(1, &texname);	
            glBindTexture(GL_TEXTURE_3D, texname);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
            glTexImage3D(GL_TEXTURE_3D, 0, GL_INTENSITY8, width, height, depth, 0, GL_RED, GL_UNSIGNED_BYTE, &texels[0]);
        }
        
        return texname;
    }
    
    
    int WireframeRenderer::maximum_face_valency(const Manifold& m)
    {
        int max_val = 0;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
            max_val = max(max_val, no_edges(m, *f));
        return max_val;
    }
    
    WireframeRenderer::WireframeRenderer(HMesh::Manifold& m, bool smooth): idbuff_renderer(0)
    {
        if(GLEW_EXT_geometry_shader4 && maximum_face_valency(m) > 3)
        {
            GLint viewp[4];
            glGetIntegerv(GL_VIEWPORT,viewp);
            idbuff_renderer = new IDBufferWireframeRenderer(viewp[2], viewp[3], m);
        }
        else
        {
            glNewList(display_list,GL_COMPILE);
            if(GLEW_EXT_geometry_shader4)
                draw_triangles_in_wireframe(m,smooth, Vec3f(1,0,0));				
            else
                draw_wireframe_oldfashioned(m,smooth, Vec3f(1,0,0));
            glEndList();
        }
    }
    
    void WireframeRenderer::draw()
    {
        if(idbuff_renderer)
        {
            glEnable(GL_LIGHTING);
            idbuff_renderer->draw(Vec3f(1,0,0),Vec3f(1));
            glDisable(GL_LIGHTING);
        }
        else
            glCallList(display_list);
    }
    
    void SimpleShaderRenderer::init_shaders(const std::string& vss, 
                                            const std::string& fss)
    {
        vs = create_glsl_shader(GL_VERTEX_SHADER, vss);
        print_glsl_program_log(vs);
        
        fs = create_glsl_shader(GL_FRAGMENT_SHADER, fss);
        print_glsl_program_log(fs);
        
        prog = glCreateProgram();
        
        if(vs) glAttachShader(prog, vs);
        if(fs) glAttachShader(prog, fs);
        
        glLinkProgram(prog);
        print_glsl_program_log(prog);
        
    }
    
    void SimpleShaderRenderer::compile_display_list(const Manifold& m, bool smooth)
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glNewList(display_list,GL_COMPILE);
        GLGraphics::draw(m, smooth);
        glEndList();	
        glUseProgram(old_prog);
    }
    
    void SimpleShaderRenderer::draw()
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glCallList(display_list);
        glUseProgram(old_prog);
    }
    
    const string NormalRenderer::vss =
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string NormalRenderer::fss =
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 l = normalize(-v);\n"
    "	vec3 e = l;\n"
    "	vec3 r = normalize(2.0*dot(l, n)*n - l);\n"
    "	\n"
    "	vec4 a = vec4(0.0,0.1,.3,1.0);\n"
    "   float dot_ln = abs(dot(l, n));\n"
    "	vec4 d = vec4(0.7)*dot_ln;\n"
    "	vec4 s = vec4(0.3)*smoothstep(0.98,0.9999,dot(r, e));\n"
    "	\n"
    "	gl_FragColor = d+s;\n"
    "}\n";
    
    const string DebugRenderer::vss =
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 c;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "   c = gl_Color.rgb;\n"
    "}\n";
    
    const string DebugRenderer::fss =
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 c;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 l = normalize(-v);\n"
    "	vec3 e = l;\n"
    "	vec3 r = normalize(2.0*dot(l, n)*n - l);\n"
    "	\n"
    "	vec4 a = vec4(0.0,0.1,.3,1.0);\n"
    "   float dot_ln = abs(dot(l, n));\n"
    "	vec4 d = vec4(c,1) * 0.7 * dot_ln;\n"
    "	vec4 s = vec4(c,1) * 0.3 * smoothstep(0.98,0.9999,dot(r, e));\n"
    "	\n"
    "	gl_FragColor =  vec4(c,1);\n"
    "}\n";
    
    HMesh::VertexAttributeVector<CGLA::Vec3f> DebugRenderer::vertex_colors;
    HMesh::HalfEdgeAttributeVector<CGLA::Vec3f> DebugRenderer::edge_colors;
    HMesh::FaceAttributeVector<CGLA::Vec3f> DebugRenderer::face_colors;

    
    void DebugRenderer::compile_display_list(const HMesh::Manifold& m, bool smooth, float rad)
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glNewList(display_list,GL_COMPILE);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1,1);
        for(FaceID f: m.faces()){
            Vec3f c = face_colors[f];
            glColor3f(c[0], c[1], c[2]);
            if(!smooth)
                glNormal3dv(normal(m, f).get());
            if(no_edges(m, f)== 3)
                glBegin(GL_TRIANGLES);
            else
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n = normal(m, w.vertex());
                if(smooth)
                    glNormal3dv(n.get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glLineWidth(2);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glBegin(GL_LINES);
        for(auto hid: m.halfedges())
        {
            Walker w = m.walker(hid);
            Vec3f c = edge_colors[hid];
            glColor3fv(c.get());
            glNormal3dv(normal(m, w.opp().vertex()).get());
            glVertex3dv(m.pos(w.opp().vertex()).get());
            glNormal3dv(normal(m, w.vertex()).get());
            glVertex3dv(m.pos(w.vertex()).get());
        }
        glEnd();
        glLineWidth(1);
        Vec3d c;
        float r;
        bsphere(m, c, r);
        r *= rad;
        for(auto vid : m.vertices())
        {
            Vec3d p = m.pos(vid);
            Vec3f c = vertex_colors[vid];
            glColor3f(c[0], c[1], c[2]);
            glPushMatrix();
            glTranslated(p[0], p[1], p[2]);
            glScalef(r, r, r);
            draw_ball();
            glPopMatrix();
        }
        glEnd();
        glEndList();
        glUseProgram(old_prog);
    }
    
    const string ReflectionLineRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    
    const string ReflectionLineRenderer::fss = 
    "uniform float detail;\n"
    "\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	// calculate the reflection\n"
    "	vec3 r = normalize(2.0*dot(-v, n)*n + v);\n"
    "	vec3 viewer_lightdir = vec3(0, 0, 1.0);\n"
    "   float diff  = dot(n,viewer_lightdir);\n"
    "	\n"
    "	vec2 r2 = normalize(vec2(r[0], r[2]));\n"
    "	vec2 x = vec2(1, 0);\n"
    "	float angle = acos(dot(r2, x));\n"
    "	\n"
    "	// decide if we hit a white or black ring, based on y value\n"
    "	gl_FragColor = diff * vec4(1.0) + smoothstep(0.9, 0.95,cos(13.0*angle)) * vec4(-1.0);\n"
    "}\n";
    
    const string IsophoteLineRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    
    const string IsophoteLineRenderer::fss = 
    "uniform float detail;\n"
    "\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 viewer_lightdir = vec3(0, 0, 1.0);\n"
    "	vec3 isophote_lightdir = viewer_lightdir;\n"
    "	float angle = acos(dot(n, isophote_lightdir));\n"
    "   float diff  = dot(n,viewer_lightdir);\n"
    "	\n"
    "	// decide if we hit a white or black ring, based on y value\n"
    "	gl_FragColor = diff * vec4(1.0) + smoothstep(0.9,0.95,cos(20.0*angle)) * vec4(-1.0);\n"
    "}\n";
    
    const string ToonRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string ToonRenderer::fss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 l = normalize(-v);\n"
    "	vec3 e = l;\n"
    "	vec3 r = normalize(2.0*dot(l, n)*n - l);\n"
    "	\n"
    "	vec4 a = vec4(0.0,0.1,.3,1.0);\n"
    "   float dot_ln = abs(dot(l, n));\n"
    "	vec4 d = vec4(0.7,0.7,0.0,1.0) * 0.25 * (smoothstep(0.23,0.25,dot_ln)+smoothstep(0.45,0.47,dot_ln)+smoothstep(0.7,0.72,dot_ln)+smoothstep(0.9,0.92,dot_ln));\n"
    "	vec4 s = vec4(0.5,0.3,0.4,1.0)*smoothstep(0.96,0.98,dot(r, e));\n"
    "	\n"
    "	gl_FragColor =  d+s;\n"
    "}\n";
    
    
    void GlazedRenderer::compile_display_list(const HMesh::Manifold& m, bool smooth)
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glNewList(display_list,GL_COMPILE);
        glBindTexture(GL_TEXTURE_3D, get_noise_texture_id());
        glUniform1iARB(glGetUniformLocationARB(prog, "noise_tex"),0);
        float r;
        Vec3d c;
        bsphere(m, c, r);
        glUniform1fARB(glGetUniformLocationARB(prog, "noise_scale"),12.0/r);
        GLGraphics::draw(m, smooth);
        glEndList();
        glUseProgram(old_prog);

    }

    
    const string GlazedRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "   v_obj = gl_Vertex.xyz;\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n"
    "\n";
    
    const string GlazedRenderer::fss =
    "uniform sampler3D noise_tex;\n"
    "uniform float noise_scale;\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "vec4 glazed_shader(vec4 mat_col,  vec4 light_col, vec3 light_dir)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 e = normalize(-v);\n"
    "	vec3 r = normalize(2.0*dot(e, n)*n - e);\n"
    "	float d = max(0.05,dot(light_dir, n));\n"
    "	vec4 diff = mat_col * light_col *d; 	\n"
    "	vec4 refl = smoothstep(0.7,0.75,dot(r,light_dir)) * light_col;\n"
    "	return 0.15*refl + diff;\n"
    "}\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	vec4 mat_col = vec4(0.9,1.0,0.4,1.0) +  vec4(-0.1,-0.1,0.12,0.0) * texture3D(noise_tex, noise_scale*v_obj).x\n"
    " + vec4(0.05) * texture3D(noise_tex, 500.0*v_obj).x;\n"
    "	\n"
    "	vec3 light0_dir = vec3(0.0,1.0,0.0);\n"
    "	vec4 light0_col = vec4(0.7,0.9,1.0,1.0);\n"
    "	\n"
    "	vec3 light1_dir = vec3(0.0,0.0,1.0);\n"
    "	vec4 light1_col = vec4(1.0,1.0,0.7,1.0);\n"
    "	\n"
    "	gl_FragColor = \n"
    "	0.5*glazed_shader(mat_col, light0_col, light0_dir)+\n"
    "	0.5*glazed_shader(mat_col, light1_col, light1_dir);\n"
    "	\n"
    "	gl_FragColor.a = 1.0;\n"
    "}\n";
    
    
    const string ScalarFieldRenderer::vss =
    "	attribute float scalar;\n"
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		s=scalar;\n"
    "	}\n";
    
    const string ScalarFieldRenderer::fss =
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	uniform float scalar_min;\n"
    "	uniform float scalar_max;\n"
    "   uniform float gamma;\n"
    "   uniform int use_shading;\n"
    "   uniform int use_stripes;\n"
    "   uniform int color_signed;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    " vec4 rainbow(float f) {\n"
    "    const float dx = 0.8;\n"
    "    float g = (6.-2.*dx)*f+dx;\n"
    "    float R = max(0.0,(3.-abs(g-4.)-abs(g-5.))/2.0);\n"
    "    float G = max(0.0,(4.-abs(g-2.)-abs(g-4.))/2.0);\n"
    "    float B = max(0.0,(3.-abs(g-1.)-abs(g-2.))/2.0);\n"
    "    return vec4(R,G,B,0.0);\n"
    " }"
    " float stripe(float x) {return fract(x*20.0)>0.5?1.0:0.0;}"
    " float saw(float x) {const float y=0.05; return  (x - y * floor(x/y))/y;}"
    "	void main()\n"
    "	{\n"
    "       vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		\n"
    "       if(color_signed==1) {\n"
    "		  float s_norm = s/max(abs(scalar_max),abs(scalar_min));\n"
    "		  gl_FragColor = s_norm*vec4(1,0,-1,0);\n"
    "         if(use_stripes==1) gl_FragColor += vec4(.01,.01,.01,0)*stripe(s_norm);"
    "       } else {\n"
    "		  float s_norm = (s-scalar_min)/(scalar_max-scalar_min);\n"
    "		  gl_FragColor = rainbow(s_norm);\n"
    "         if(use_stripes==1) gl_FragColor -= vec4(.1,.1,.1,0)*stripe(s_norm);"
    "       }\n"
    "       if(use_shading==1) gl_FragColor *= 0.5+0.5*dot_ln;\n"
    "       gl_FragColor.r = pow(max(0.0,gl_FragColor.r), 1.0/gamma);\n"
    "       gl_FragColor.g = pow(max(0.0,gl_FragColor.g), 1.0/gamma);\n"
    "       gl_FragColor.b = pow(max(0.0,gl_FragColor.b), 1.0/gamma);\n"
    "	}\n";
    
    void ScalarFieldRenderer::compile_display_list(const HMesh::Manifold& m, bool smooth,
                                        HMesh::VertexAttributeVector<double>& field,
                                        double min_val, double max_val,
                                        float gamma,
                                        int use_stripes,
                                        int color_signed,
                                        int use_shading)
    {
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint scalar_attrib = glGetAttribLocation(prog, "scalar");
        glUniform1fARB(glGetUniformLocationARB(prog, "scalar_max"), max_val);
        glUniform1fARB(glGetUniformLocationARB(prog, "scalar_min"), min_val);
        glUniform1iARB(glGetUniformLocationARB(prog, "use_shading"), use_shading);
        glUniform1iARB(glGetUniformLocationARB(prog, "use_stripes"), use_stripes);
        glUniform1iARB(glGetUniformLocationARB(prog, "color_signed"), color_signed);
        
        //    static float& gamma = CreateCVar("display.scalar_field_renderer.gamma",2.2f);
        glUniform1fARB(glGetUniformLocationARB(prog, "gamma"), gamma);
        glNewList(display_list,GL_COMPILE);
        
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(!smooth)
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3)
                glBegin(GL_TRIANGLES);
            else
                glBegin(GL_POLYGON);
            
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                if(smooth)
                    glNormal3dv(n.get());
                glVertexAttrib1d(scalar_attrib, field[w.vertex()]);
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();
        glUseProgram(old_prog);
        
    }
    
    HMesh::VertexAttributeVector<CGLA::Vec2f> CheckerBoardRenderer::param;

    const string CheckerBoardRenderer::vss =
    "	attribute vec2 param;\n"
    "	varying vec3 _normal;\n"
    "	varying vec2 uv;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		uv=param;\n"
    "	}\n";
    
    const string CheckerBoardRenderer::fss =
    "	varying vec3 _normal;\n"
    "	varying vec2 uv;\n"
    "   const float pi = 3.14159265359;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    "	void main()\n"
    "	{\n"
    "       vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		vec2 rt = uv;//vec2(length(uv),atan(uv.y, uv.x));\n"
    "		float stripe_signal = smoothstep(-0.001,0.001,sin(2.0*pi*rt.x)*sin(10.0*rt.y));\n"
    "		\n"
    "		gl_FragColor = dot_ln * vec4(0.35,0.25,0.5,0);\n"
   "		gl_FragColor.rgb += 0.7*stripe_signal;\n"
    "	}\n";
    
    void CheckerBoardRenderer::compile_display_list(const HMesh::Manifold& m, bool smooth)
    {
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint param_attrib = glGetAttribLocation(prog, "param");
        glNewList(display_list,GL_COMPILE);

        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(!smooth)
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3)
                glBegin(GL_TRIANGLES);
            else
                glBegin(GL_POLYGON);
            
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                if(smooth)
                    glNormal3dv(n.get());
                glVertexAttrib2fv(param_attrib, param[w.vertex()].get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();
        glUseProgram(old_prog);
        
    }



    const string AmbientOcclusionRenderer::vss =
    "	attribute float scalar;\n"
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		s=scalar;\n"
    "	}\n";
    
    const string AmbientOcclusionRenderer::fss = 	
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	uniform float scalar_max;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    "	void main()\n"
    "	{\n"
    "   vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		\n"
    "		float s_norm = min(1.0,s/scalar_max+1.0);\n"
    "		\n"
    "		gl_FragColor = s_norm * vec4(1.0);\n"
    "       gl_FragColor *= dot_ln;\n"
    "       gl_FragColor.r = pow(gl_FragColor.r, 1.0);\n"
    "       gl_FragColor.g = pow(gl_FragColor.g, 1.0);\n"
    "       gl_FragColor.b = pow(gl_FragColor.b, 1.0);\n"
    "	}\n";
    
    void AmbientOcclusionRenderer::compile_display_list(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& field, double max_val)
    {	
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint scalar_attrib = glGetAttribLocation(prog, "scalar");
        glUniform1fARB(glGetUniformLocationARB(prog, "scalar_max"), max_val);
        
        glNewList(display_list,GL_COMPILE);
        
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){

            if(no_edges(m, *f)== 3)
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
            {
                Vec3d n(normal(m, w.vertex()));
                glNormal3dv(n.get());
                glVertexAttrib1d(scalar_attrib, field[w.vertex()]);
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();	
        glUseProgram(old_prog);
        
    }
    
    
    void LineFieldRenderer::compile_display_list(const HMesh::Manifold& m,HMesh::VertexAttributeVector<CGLA::Vec3d>& lines)
    {
        float r;
        Vec3d c;
        bsphere(m, c, r);
        float noise_scale = 10.0f/r;
        float line_scale = 0.02f;
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glNewList(display_list,GL_COMPILE);
        glUniform1fARB(glGetUniformLocationARB(prog, "line_scale"),line_scale);
        glUniform1fARB(glGetUniformLocationARB(prog, "noise_scale"),noise_scale);
        glUniform1iARB(glGetUniformLocationARB(prog, "noise_tex"),0);
        GLuint direction = glGetAttribLocation(prog, "direction");	
        glBindTexture(GL_TEXTURE_3D, get_noise_texture_id());


        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(no_edges(m, *f) == 3)
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            Vec3d n(normal(m, *f));
            Vec3d d0 = lines[m.walker(*f).vertex()];
            d0 = normalize(d0-n*dot(n,d0));
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                glNormal3dv(n.get());
                
                Vec3d d = lines[w.vertex()];
                d = normalize(d-n*dot(n,d));
                if(dot(d,d0)<0) d=-d;
                glVertexAttrib3dv(direction, d.get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }

        glBindTexture(GL_TEXTURE_3D, 0);
        glEndList();	
        glUseProgram(old_prog);
        
    }
    
    
    const string LineFieldRenderer::vss = 
    "attribute vec3 direction;\n"
    "varying vec3 _n;\n"
    "varying vec3 dir_obj;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "   v_obj = gl_Vertex.xyz;\n"
    "	dir_obj = direction;\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string LineFieldRenderer::fss =
    "uniform sampler3D noise_tex;\n"
    "uniform float line_scale;\n"
    "uniform float noise_scale;\n"
    "varying vec3 _n;\n"
    "varying vec3 dir_obj;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "float tex(vec3 p) {return smoothstep(0.2,0.4,texture3D(noise_tex, p).x);}\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "   vec3 d = normalize(dir_obj);\n"
    "   float I = "
    "             tex(noise_scale*v_obj + 6.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 6.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj + 5.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 5.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj + 4.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 4.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj + 3.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 3.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj + 2.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 2.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj + 1.0*line_scale*d) + \n"
    "             tex(noise_scale*v_obj - 1.0*line_scale*d) + \n"
    "			  tex(noise_scale*v_obj); \n"
    "	\n"
    "   float diff = max(0.0,dot(n,vec3(0.0, 0.0, 1.0)));\n"
    "	gl_FragColor.rgb = vec3(diff*I/13.0);\n"
    "	gl_FragColor.a = 1.0;\n"
    "}\n";
    
    GLuint HarmonicsRenderer::prog_P0 = 0;
    GLGraphics::Console::variable<float> HarmonicsRenderer::display_harmonics_time;
    GLGraphics::Console::variable<int> HarmonicsRenderer::display_harmonics_diffuse;
    GLGraphics::Console::variable<int> HarmonicsRenderer::display_harmonics_highlight;

    GLGraphics::Console::variable<int> HarmonicsRenderer::display_harmonics_e0;
    GLGraphics::Console::variable<int> HarmonicsRenderer::display_harmonics_e1;
    
    string vss =
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"	\n"
	"	\n"
	"	attribute float eigenvalue;\n"
	"	attribute float eigenvalue2;\n"
	"	varying vec3 normal;\n"
	"	varying float eig;\n"
	"	varying float eig2;\n"
	"	\n"
	"	void main(void)\n"
	"	{\n"
	"		gl_Position =  ftransform();\n"
	"		normal = normalize(gl_NormalMatrix * gl_Normal);\n"
	"		eig = eigenvalue;\n"
	"		eig2 = eigenvalue2;\n"
	"	}\n";
	
	string fss =
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"	\n"
	"	varying vec3 normal;\n"
	"	varying float eig;\n"
	"	varying float eig2;\n"
	"	uniform float eig_max;\n"
	"	uniform float eig_max2;\n"
	"	uniform bool do_highlight;\n"
	"	uniform bool do_diffuse;\n"
	"	const vec3 light_dir = vec3(0,0,1);\n"
	"	\n"
	" float basef(float x) {return max(0.0,min(1.0,2.0-4.0*abs(x)));\n}"
	"	void main()\n"
	"	{\n"
	"		float dot_ln = max(0.0,dot(light_dir, normal));\n"
	"		\n"
	"		float eig_norm = eig/eig_max;\n"
	"		float stripe_signal = 250 * eig_norm;\n"
	//"		vec4 stripe_col = abs(stripe_signal) < 3.14 ? vec4(1,1,0,0) : vec4(.4,.4,.4,0);\n"
	"		vec4 stripe_col = -vec4(.4,.4,.4,0);\n"
	"		\n"
	"       float alpha = (1.0-eig_norm) * 2.0 * 3.1415926;\n"
	"       float offs = 2.0*3.1415/3.0;\n"
	"		gl_FragColor = vec4(0,0,1,0)*basef(eig_norm)+vec4(0,1,0,0)*basef(eig_norm-0.5)+vec4(1,0,0,0)* basef(eig_norm-1.0);\n"
	"		if(do_diffuse)   gl_FragColor *= dot_ln;\n"
	"		if(do_highlight) gl_FragColor += dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*vec4(.5,.5,.5,0);\n"
	"		gl_FragColor -= vec4(.4,.4,.4,.4)*smoothstep(0.2,0.6,cos(stripe_signal));\n"
	"	}\n";
    

    HarmonicsRenderer::HarmonicsRenderer(HMesh::Manifold& _m, HMesh::Harmonics* _h, GLGraphics::Console& cs): m(&_m), h(_h)
    {
        if (prog_P0 == 0) {
            string shader_path = "/Users/jab/GEL/apps/MeshEdit/";
            GLuint vs = create_glsl_shader(GL_VERTEX_SHADER, vss);
            GLuint fs = create_glsl_shader(GL_FRAGMENT_SHADER, fss);
            
            // Create the program
            prog_P0 = glCreateProgram();
            
            // Attach all shaders
            if(vs) glAttachShader(prog_P0, vs);
            if(fs) glAttachShader(prog_P0, fs);
            
            // Link the program object and print out the info log
            glLinkProgram(prog_P0);
            print_glsl_program_log(prog_P0);
            
            // Install program object as part of current state
            glUseProgram(0);
            
            
            display_harmonics_diffuse.reg(cs, "display.harmonics.diffuse", "");
            display_harmonics_time.reg(cs, "display.harmonics.time", "");
            display_harmonics_highlight.reg(cs, "display.harmonics.highlight", "");
            display_harmonics_e0.reg(cs,"display.harmonics.e0","");
            display_harmonics_e1.reg(cs,"display.harmonics.e1","");
        }
        draw_esum();
    }
    
    
    void HarmonicsRenderer::parse_key(unsigned char key)
    {
        switch(key) {
            case '+':
                display_harmonics_time = display_harmonics_time+0.001;
                break;
            case '-':
                display_harmonics_time = display_harmonics_time-0.001;
                break;
            case 'd':
                display_harmonics_diffuse = !display_harmonics_diffuse;
                break;
            case 'h':
                display_harmonics_highlight = !display_harmonics_highlight;
                break;
        }
        
    }
    
    
    
    
    void HarmonicsRenderer::draw_adf()
    {
        VertexAttributeVector<double> F;
        double F_max = h->compute_adf(F, display_harmonics_time);
        cout << "F max" <<  F_max << endl;
        
        glNewList(display_list, GL_COMPILE);
        glUseProgram(prog_P0);
        glUniform1f(glGetUniformLocation(prog_P0,"eig_max"),F_max);//2*M_PI);
        glUniform1i(glGetUniformLocation(prog_P0,"do_diffuse"),display_harmonics_diffuse);
        glUniform1i(glGetUniformLocation(prog_P0,"do_highlight"),display_harmonics_highlight);
        GLuint attrib = glGetAttribLocationARB(prog_P0, "eigenvalue");
        
        glFrontFace(GL_CW);
        for(FaceIDIterator f = m->faces_begin(); f != m->faces_end(); ++f){
            glBegin(GL_TRIANGLES);
            for(Walker w = m->walker(*f); !w.full_circle(); w = w.circulate_face_cw()){
                glVertexAttrib1f(attrib,F[w.vertex()]);
                glNormal3dv(normal(*m, w.vertex()).get());
                glVertex3dv(m->pos(w.vertex()).get());
            }
            glEnd();
        }
        glFrontFace(GL_CCW);
        glUseProgram(0);
        glEndList();
    }

    void HarmonicsRenderer::draw_esum()
    {
        VertexAttributeVector<double> F;
        double F_max = h->compute_esum(F, display_harmonics_e0, display_harmonics_e1);
        cout << "F max" <<  F_max << endl;
        
        glNewList(display_list, GL_COMPILE);
        glUseProgram(prog_P0);
        glUniform1f(glGetUniformLocation(prog_P0,"eig_max"),F_max);//2*M_PI);
        glUniform1i(glGetUniformLocation(prog_P0,"do_diffuse"),display_harmonics_diffuse);
        glUniform1i(glGetUniformLocation(prog_P0,"do_highlight"),display_harmonics_highlight);
        GLuint attrib = glGetAttribLocationARB(prog_P0, "eigenvalue");
        
        glFrontFace(GL_CW);
        for(FaceIDIterator f = m->faces_begin(); f != m->faces_end(); ++f){
            glBegin(GL_TRIANGLES);
            for(Walker w = m->walker(*f); !w.full_circle(); w = w.circulate_face_cw()){
                glVertexAttrib1f(attrib,F[w.vertex()]);
                glNormal3dv(normal(*m, w.vertex()).get());
                glVertex3dv(m->pos(w.vertex()).get());
            }
            glEnd();
        }
        glFrontFace(GL_CCW);
        glUseProgram(0);
        glEndList();
    }


    
}


