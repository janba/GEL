/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "IDBufferWireFrameRenderer.h"

#include "../CGLA/Vec4f.h"
#include "../CGLA/Vec2f.h"
#include "../CGLA/Vec3f.h"
#include "../HMesh/Manifold.h"
#include "../HMesh/AttributeVector.h"

#include "../GLGraphics/draw.h"

#include "glsl_shader.h"

using namespace std;
using namespace CGLA;
using namespace GLGraphics;
using namespace HMesh;

namespace
{
    string line_atten_frag = 
        "uniform int ATTEN_MODE;\n"
        "uniform float THICKNESS;\n"
        "uniform sampler2DRect IDMAP;\n"
        "uniform float TRANSITION;\n"
        "\n"
        "varying vec3 _id;\n"
        "varying vec2 p_2d_0;\n"
        "varying vec2 p_2d_1;\n"
        "\n"
        "const int LINEAR_ATTENUATION = 1;\n"
        "const int HERMITE_ATTENUATION = 2;\n"
        "\n"
        "const float HERMITE_COL_ATTEN_BIAS = 0.25;\n"
        "const float LINEAR_COL_ATTEN_BIAS = 0.4;\n"
        "const float HERMITE_ATTEN_BEGIN = 0.75;\n"
        "const float HERMITE_ATTEN_END = 1.0;\n"
        "\n"
        "void main(void)\n"
        "{\n"
        "	vec2 dir = (p_2d_1-p_2d_0);\n"
        "	float sq_len01 = dot(dir,dir);\n"
        "	vec2 v0 = (gl_FragCoord.xy)-p_2d_0;\n"
        "	vec2 v1 = (gl_FragCoord.xy)-p_2d_1;\n"
        "	float t = dot(dir, v0);\n"
        "	float d;\n"
        "	\n"
        "	if(t<=0.0)\n"
        "			d = length(v0);\n"
        "	else if(t>= sq_len01)\n"
        "			d = length(v1);\n"
        "	else\n"
        "			d = length(dir * t / sq_len01 - v0);\n"
        "	\n"
        "	vec3 iddiff = texture2DRect(IDMAP,gl_FragCoord.xy).xyz - _id;\n"
        " 	if(dot(iddiff, iddiff)<1e-6)\n"
        "		gl_FragDepth = 0.0;\n"
        "	else	\n"
        "		gl_FragDepth = gl_FragCoord.z;\n"
        "	\n"
        "	\n"
        "	float t1 = THICKNESS;\n"
        "	if(ATTEN_MODE == HERMITE_ATTENUATION)\n"
        "		{\n"
        "			float width_atten = smoothstep(HERMITE_ATTEN_END,HERMITE_ATTEN_BEGIN,\n"
        "																		 gl_FragCoord.z);\n"
        "			t1 *= width_atten;\n"
        "		}\n"
        "	else if(ATTEN_MODE == LINEAR_ATTENUATION)\n"
        "		{\n"
        "			float width_atten = 1.0-gl_FragCoord.z;\n"
        "			t1 *= width_atten;\n"
        "		}			\n"
        "	float t2 = t1+TRANSITION;\n"
        "\n"
        "	float I;\n"
        "	if(d<t1) \n"
        "			I = 1.0;\n"
        "	else\n"
        "	{\n"
        "			float x = (d-t1);\n"
        "			I = exp2(-x*x*2);\n"
        "	}\n"
        " 	gl_FragColor.rgb = gl_Color.rgb;\n"
        "	gl_FragColor.a = I;\n"
        "}\n";

    string line_frag = 
        "uniform sampler2DRect IDMAP;\n"
        "uniform float TRANSITION;\n"
        "	\n"
        "varying vec3 _id;\n"
        "varying vec2 p_2d_0;\n"
        "varying vec2 p_2d_1;\n"
        "\n"
        "void main(void)\n"
        "{\n"
        "	vec2 dir = (p_2d_1-p_2d_0);\n"
        "	float sq_len01 = dot(dir,dir);\n"
        "	vec2 v0 = gl_FragCoord.xy-p_2d_0;\n"
        "	vec2 v1 = gl_FragCoord.xy-p_2d_1;\n"
        "	float t = dot(dir, v0);\n"
        "	float sq_d;\n"
        "	\n"
        "	if(t<=0.0)\n"
        "		sq_d = dot(v0,v0);\n"
        "	else if(t>= sq_len01)\n"
        "		sq_d = dot(v1,v1);\n"
        "	else\n"
        "		{\n"
        "			float area = (v0.x*v1.y - v0.y * v1.x);\n"
        "			sq_d = area*area/sq_len01;\n"
        "		}\n"
        "\n"
        "	vec3 iddiff = texture2DRect(IDMAP,gl_FragCoord.xy).xyz -_id;\n"
        "	if(dot(iddiff, iddiff)<1e-6)\n"
        "		gl_FragDepth = 0.0;\n"
        " 	else	\n"
        " 		gl_FragDepth = gl_FragCoord.z;\n"
        "	\n"
        "	gl_FragColor.a = exp2(-sq_d*2.0);\n"
        "	gl_FragColor.rgb = vec3(1,0,0);\n"
        "}\n";

    string line_vert = 
        "uniform float THICKNESS;\n"
        "uniform vec2 WIN_SCALE;\n"
        "uniform vec2 INV_WIN_SCALE;\n"
        "uniform float TRANSITION;\n"
        "\n"
        "attribute vec4 id;\n"
        "attribute vec4 opp_vertex;\n"
        "attribute vec2 displace;\n"
        "\n"
        "varying vec3 _id;\n"
        "varying vec2 p_2d_0;\n"
        "varying vec2 p_2d_1;\n"
        "\n"
        "\n"
        "void main(void)\n"
        "{\n"
        "	vec4 p0 = gl_ModelViewProjectionMatrix * opp_vertex;\n"
        "	float w0 = p0.w; \n"
        "\n"
        "	vec4 p1 = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
        "	float w1 = p1.w;\n"
        "		\n"
        "	p_2d_0 = (p0.xy/w0+vec2(1,1)) * WIN_SCALE;\n"
        "	p_2d_1 = (p1.xy/w1+vec2(1,1)) * WIN_SCALE;\n"
        "	\n"
        "	vec2 a = normalize(p_2d_1 - p_2d_0);\n"
        "	vec2 a_h = vec2(-a.y, a.x);\n"
        "	vec2 scale = (TRANSITION+THICKNESS)*INV_WIN_SCALE;\n"
        "\n"
        "	gl_Position = (displace.x>0.0) ? p1 : p0;\n"
        "	vec2 offset = a * displace.x + a_h * displace.y;\n"
        "	gl_Position.xy += scale * gl_Position.w * offset;\n"
        "\n"
        "	_id = id.rgb;\n"
        "}\n";

}

namespace GLGraphics
{
    IDBufferWireframeRenderer::~IDBufferWireframeRenderer()
    {
        glDeleteShader(vs);	
        glDeleteShader(fs);
        glDeleteProgram(line_prog);
        glDeleteTextures(1, &idmap);
        glDeleteBuffers(1, &vertex_buffername);
        glDeleteBuffers(1, &colors_buffername);

        glDeleteBuffers(1, &line_id_attrib);
        glDeleteBuffers(1, &line_vertex_pos);
        glDeleteBuffers(1, &line_disp_attrib);			
        glDeleteBuffers(1, &line_opp_attrib);
    }

    IDBufferWireframeRenderer::IDBufferWireframeRenderer(int _XSZ, int _YSZ,
        HMesh::Manifold& _mesh, 
        float _thickness, 
        float _transition, 
        int atten_mode): 
    mesh(&_mesh), XSZ(_XSZ), YSZ(_YSZ), thickness(_thickness), transition(_transition)
    {

        if(atten_mode == 0 && thickness == 0.0)
        {
            vs = create_glsl_shader(GL_VERTEX_SHADER, line_vert);
            fs = create_glsl_shader(GL_FRAGMENT_SHADER, line_frag);
        }
        else
        {
            vs = create_glsl_shader(GL_VERTEX_SHADER, line_vert);
            fs = create_glsl_shader(GL_FRAGMENT_SHADER, line_atten_frag);
        }

        line_prog = glCreateProgram();
        glAttachShader(line_prog, vs);
        glAttachShader(line_prog, fs);
        glLinkProgram(line_prog);
        glUseProgram(line_prog);

        glUniform1f(glGetUniformLocation(line_prog,"THICKNESS"), 
            thickness);
        glUniform1f(glGetUniformLocation(line_prog,"TRANSITION"), 
            transition);
        glUniform1i(glGetUniformLocation(line_prog,"ATTEN_MODE"),
            atten_mode);
        glUniform2f(glGetUniformLocation(line_prog,"WIN_SCALE"), 
            XSZ/2.0, YSZ/2.0);
        glUniform2f(glGetUniformLocation(line_prog,"INV_WIN_SCALE"), 
            2.0/XSZ, 2.0/YSZ);
        glUniform1i(glGetUniformLocation(line_prog,"IDMAP"), 0);

        id_attrib = glGetAttribLocation(line_prog, "id");
        popp_attrib = glGetAttribLocation(line_prog, "opp_vertex");
        disp_attrib = glGetAttribLocation(line_prog, "displace");
        glUseProgram(0);

        glGenTextures(1, &idmap);
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, idmap);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,GL_TEXTURE_WRAP_S,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,GL_TEXTURE_WRAP_T,GL_CLAMP);

        //create the texture
        glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA8, XSZ, YSZ,
            0, GL_RGB, GL_FLOAT, 0);


        glGenBuffers(1, &vertex_buffername);
        glGenBuffers(1, &colors_buffername);

        glGenBuffers(1, &line_id_attrib);
        glGenBuffers(1, &line_vertex_pos);
        glGenBuffers(1, &line_disp_attrib);
        glGenBuffers(1, &line_opp_attrib);

        triangles = static_cast<int>(mesh->no_faces());
        vector<Vec3f> verts;
        vector<Vec3f> cols;

        unsigned int k = 0;
        for(FaceIDIterator f = mesh->faces_begin(); f != mesh->faces_end(); ++f, ++k){
            Vec3uc idv(id_get(k));
            Vec3f idvec(idv[0]/255.0, idv[1]/255.0, idv[2]/255.0);
            for(Walker w = mesh->walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                cols.push_back(idvec);
                verts.push_back(Vec3f(mesh->pos(w.vertex())));
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffername);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*verts.size(),
            (float*)&verts[0],GL_STATIC_DRAW);



        glBindBuffer(GL_ARRAY_BUFFER, colors_buffername);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*cols.size(),
            (float*)&cols[0],GL_STATIC_DRAW);


        vector<Vec3f> line_ids;
        vector<Vec3f> vertex_positions;
        vector<Vec2f> displacements;
        vector<Vec3f> opposite_positions;


        quads = 0;
        unsigned int i = 0;
        for(FaceIDIterator f = mesh->faces_begin(); f != mesh->faces_end(); ++f,++i){
            for(Walker w = mesh->walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                ++quads;
                Vec3uc idv(id_get(i));
                Vec3f v0(mesh->pos(w.next().vertex()));
                Vec3f v1(mesh->pos(w.next().opp().vertex()));
                Vec3f idvec(idv[0]/255.0, idv[1]/255.0, idv[2]/255.0);

                line_ids.push_back(idvec);
                opposite_positions.push_back(v0);
                displacements.push_back(Vec2f(1,-1));
                vertex_positions.push_back(v1);


                line_ids.push_back(idvec);
                opposite_positions.push_back(v0);
                displacements.push_back(Vec2f(1, 1));
                vertex_positions.push_back(v1);


                line_ids.push_back(idvec);
                opposite_positions.push_back(v0);
                displacements.push_back(Vec2f(-1,1));
                vertex_positions.push_back(v1);


                line_ids.push_back(idvec);
                opposite_positions.push_back(v0);
                displacements.push_back(Vec2f(-1,-1));
                vertex_positions.push_back(v1);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, line_id_attrib);

        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*line_ids.size(),
            (float*)&line_ids[0],GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, line_opp_attrib);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*opposite_positions.size(),
            (float*)&opposite_positions[0],GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, line_disp_attrib);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*2*displacements.size(),
            (float*)&displacements[0],GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, line_vertex_pos);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*vertex_positions.size(),
            (float*)&vertex_positions[0],GL_STATIC_DRAW);



    }


    void IDBufferWireframeRenderer::draw(const Vec3f& color, const Vec3f& clear_color)
    {
        // push those attributes we change.
        glPushAttrib(GL_COLOR_BUFFER_BIT|
            GL_CURRENT_BIT|
            GL_TRANSFORM_BIT|
            GL_DEPTH_BUFFER_BIT);

        // Store information about whether we use lighting
        GLboolean lights_on;
        glGetBooleanv(GL_LIGHTING, &lights_on);

        // Store color information
        Vec4f current_color;
        glGetFloatv(GL_CURRENT_COLOR, &current_color[0]);

        // Store the current draw buffer 
        GLint _currentDrawbuf;
        glGetIntegerv(GL_DRAW_BUFFER, &_currentDrawbuf); 

        // Enable depth testing
        glEnable(GL_DEPTH_TEST);

        // ------------------------------
        // Pass 1: Draw the ID map
        // Each polygon has a unique ID which is coded as a colour and drawn
        // into the ID buffer
        glDisable(GL_LIGHTING);
        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffername);
        glVertexPointer(3,GL_FLOAT,0,static_cast<char*>(0));
        glEnableClientState(GL_VERTEX_ARRAY);

        glBindBuffer(GL_ARRAY_BUFFER, colors_buffername);
        glColorPointer(3,GL_FLOAT,0,static_cast<char*>(0));
        glEnableClientState(GL_COLOR_ARRAY);

        int vertex_idx = 0;
        for(FaceIDIterator f = mesh->faces_begin(); f != mesh->faces_end(); ++f){
            glBegin(GL_POLYGON);
            for(Walker w = mesh->walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
                 glArrayElement(vertex_idx++);
            glEnd();
        }
        glFinish();
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);

        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, idmap);
        glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, 0,0,0,0,XSZ, YSZ);

        // Clear color buffer but retain depth buffer.
        glClear(GL_COLOR_BUFFER_BIT);

        // Enable blending for all subsequent passes
        glEnable(GL_BLEND);

        // ------------------------------
        // PASS 3: Draw lines using ID map texture.
        // For each polygon, all edges are drawn as prefiltered
        // lines. A given fragment belonging to a line will be written
        // to the framebuffer if either of the following three conditions are met
        // - its ID matches the contents of the ID  buffer for the corresponding 
        // pixel (in the ID buffer)
        // - The ID of the corresponding pixel is 0 - meaning we are outside.
        // - The depth test is passed.
        // The final condition ensures that line pixels which are on an interior
        // contour will be drawn.
        //
        // If the line fragment is written into the framebuffer - two values are
        // actually written: The colour of the line and an alpha value which
        // corresponds to the value of the filter function.
        //
        // During this pass, blending is enabled and the blending equation is set
        // to max. Since the alpha values are the values of the filter (large if 
        // close to line) this means that we replace a pixel value with an 
        // incoming fragment if the incoming fragment is closer to a line
        // than the pixel value.
        //
        // The depth values are not changed during this pass.
        glEnable(GL_TEXTURE_RECTANGLE_ARB);
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, idmap);
        glDepthMask(GL_FALSE);
        glBlendEquation(GL_MAX);
        glUseProgram(line_prog);

        float lw;
        glGetFloatv(GL_LINE_WIDTH, &lw);
        glLineWidth(ceil(2.0*(thickness+transition)));


        glBindBuffer(GL_ARRAY_BUFFER, line_id_attrib);
        glVertexAttribPointer(id_attrib, 3,GL_FLOAT,GL_FALSE,0,static_cast<char*>(0));
        glEnableVertexAttribArray(id_attrib);


        glBindBuffer(GL_ARRAY_BUFFER, line_disp_attrib);
        glVertexAttribPointer(disp_attrib, 2,GL_FLOAT,GL_FALSE,0,static_cast<char*>(0));
        glEnableVertexAttribArray(disp_attrib);


        glBindBuffer(GL_ARRAY_BUFFER, line_opp_attrib);
        glVertexAttribPointer(popp_attrib, 3,GL_FLOAT,GL_FALSE,0,static_cast<char*>(0));
        glEnableVertexAttribArray(popp_attrib);


        glBindBuffer(GL_ARRAY_BUFFER, line_vertex_pos);
        glVertexPointer(3,GL_FLOAT,0,static_cast<char*>(0));
        glEnableClientState(GL_VERTEX_ARRAY);

        glDrawArrays(GL_QUADS, 0, quads*4);

        glDisableVertexAttribArray(id_attrib);
        glDisableVertexAttribArray(disp_attrib);
        glDisableVertexAttribArray(popp_attrib);
        glDisableClientState(GL_VERTEX_ARRAY);


        glLineWidth(lw);

        glUseProgram(0);
        glDisable(GL_TEXTURE_RECTANGLE_ARB);
        glDepthMask(GL_TRUE);

        // ------------------------------
        // Pass 4: Draw with shading
        // In this pass we draw the shaded model. At this point, the framebuffer
        // contains alpha values and line colours and also depth values.
        //
        // The depth test is set to `equal' and the shaded fragments from the 
        // filled polygons are combined with the line colours using the alpha 
        // values already stored in the frame buffer.
        //
        // The framebuffer lines along the contour also need to be blended. 
        // Hence, a screen filling quad is drawn. 
        glDepthFunc(GL_LEQUAL);
        glBlendEquation(GL_FUNC_ADD);
        glBlendFuncSeparate(GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA, 
            GL_ZERO, GL_ONE);
        if(lights_on)
        {
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
        }else
            glColor4fv(current_color.get());

        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffername);
        glVertexPointer(3,GL_FLOAT,0,static_cast<char*>(0));
        glEnableClientState(GL_VERTEX_ARRAY);

        vertex_idx = 0;
        for(FaceIDIterator f = mesh->faces_begin(); f != mesh->faces_end(); ++f){
			Vec3f n(normal(*mesh, *f));
            glNormal3fv(n.get());
            glBegin(GL_POLYGON);
            for(Walker w = mesh->walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
                 glArrayElement(vertex_idx++);
            glEnd();
        }
        glDisableClientState(GL_VERTEX_ARRAY);


        glDisable(GL_LIGHTING);


        double mvmat[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
        double prjmat[16];
        glGetDoublev(GL_PROJECTION_MATRIX, prjmat);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glColor3fv(clear_color.get());
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,1);
        glVertex3f( 1,-1,1);
        glVertex3f( 1, 1,1);
        glVertex3f(-1, 1,1);
        glEnd();



        glMatrixMode(GL_PROJECTION);
        glLoadMatrixd(prjmat);
        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixd(mvmat);

        glPopAttrib();
        if(lights_on)
            glEnable(GL_LIGHTING);

        //			cout << gluErrorString(glGetError()) << endl;
    }
}
