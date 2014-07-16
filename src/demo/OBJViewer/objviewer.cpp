// ----------------------------------------
// A simple OBJ viewer.
//
// Controls:
// - left mouse down + mouse motion : rotate
// - Scroll button and +- buttons   : zoom
// - right mouse click              : centre trackball
// - esc                            : exits
// - x,y,z buttons                  : switch trackball up axis
// - w                              : toggle wireframe on/off
// - t                              : toggle texture on/off
// - f                              : switch between vertex and face normals
// ----------------------------------------

#if (_MSC_VER >= 1200)
#pragma warning (disable: 4786)
#endif

#include <list>
#include <vector>

#include <assert.h>
#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#include <io.h>
#endif
#include <string.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <vector>

#include <GEL/Util/ArgExtracter.h>
#include <GEL/CGLA/Vec2i.h>
#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Mat4x4f.h>
#include <GEL/GLGraphics/glsl_shader.h>
#include <GEL/GLGraphics/QuatTrackBall.h>
#include <GEL/GLGraphics/draw.h>
#include <GEL/GLGraphics/SOIL.h>
#include <GEL/Geometry/TriMesh.h>
#include <GEL/Geometry/load.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/HGrid.h>

#ifdef __APPLE__
#include <GLUT/GLUT.h>
#else
#include <GL/glut.h>
#endif

using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace HMesh;
using namespace GLGraphics;

int win_size_x = 800;
int win_size_y = 800;
bool per_vertex_normals = 1;
bool redo_display_list = 1;
bool do_wireframe = false;
Vec3f line_col = Vec3f(1,0,0);
QuatTrackBall* ball;
int spin_timer = 20;
void spin(int x);
int main_window;
TriMesh mesh;
bool do_textures = true;


bool depth_pick(int x, int y,Vec3f& wp)
{
	// Enquire about the viewport dimensions
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	// Get the minimum and maximum depth values.
	float minmax_depth[2];
	glGetFloatv(GL_DEPTH_RANGE, minmax_depth);
	
	// Read a single pixel at the position of the mouse cursor.
	float depth;
	glReadPixels(x, viewport[3]-y, 1,1, GL_DEPTH_COMPONENT,
				 GL_FLOAT, (void*) &depth);
	
	// If the depth corresponds to the far plane, we clicked on the
	// background.
	if(depth == minmax_depth[1])
		return false;
	
	// The lines below copy the viewing transformation from OpenGL
	// to local variables. The call to gluLookAt must have exactly
	// the same parameters as when the scene is drawn.
	glLoadIdentity();
	ball->set_gl_modelview();
	double mvmat[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
	
	// Copy the projection matrix. We assume it is unchanged.
	double prjmat[16];
	glGetDoublev(GL_PROJECTION_MATRIX, prjmat);
	
	// Now unproject the point from screen to world coordinates.
	double ox, oy, oz;
	gluUnProject(x,viewport[3]-y,depth,
				 mvmat,prjmat,viewport,
				 &ox, &oy, &oz);
	
	wp = Vec3f(ox,oy,oz);
	
	return true;
}


void mouse_motion(int x, int y)
{
    ball->roll_ball(Vec2i(x,win_size_y-y));
}

void mouse(int btn, int state, int x, int y)
{
    y = win_size_y-y;
	if(state == GLUT_DOWN) 
	{
		if(btn == GLUT_LEFT_BUTTON) 
			ball->grab_ball(ROTATE_ACTION, Vec2i(x,y));
		else if(btn == GLUT_MIDDLE_BUTTON) 
			ball->grab_ball(ZOOM_ACTION, Vec2i(x, y));
		else if(btn == GLUT_RIGHT_BUTTON) 
			ball->grab_ball(PAN_ACTION, Vec2i(x, y));
	}
	else if(state == GLUT_UP)
		ball->release_ball();	
}

void spin(int x)
{
	ball->do_spin();
	glutTimerFunc(spin_timer, spin, 0);  
	glutPostRedisplay();
}

void setupshader()
{
	static GLuint vs,fs,prog;
	static bool was_here = false;
	if(!was_here)
	{
		was_here = true;
		const string vss = 
		"varying vec3 n;\n"
		"varying vec3 v;\n"
		"varying vec3 v_obj;\n"
		"\n"
		"void main(void)\n"
		"{\n"
		"	gl_Position = ftransform();\n"
		"   v_obj = gl_Vertex.xyz;\n"
		"	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
		"	n = normalize(gl_NormalMatrix * gl_Normal);\n"
		"}\n"
		"\n";
		
		const string fss =
		"varying vec3 n;\n"
		"varying vec3 v;\n"
		"varying vec3 v_obj;\n"
		"\n"
		"vec4 glazed_shader(vec4 mat_col,  vec4 light_col, vec3 light_dir)\n"
		"{\n"
		"	vec3 e = normalize(-v);\n"
		"	vec3 r = normalize(2.0*dot(e, n)*n - e);\n"
		"	float d = max(0.05,dot(light_dir, n));\n"
		"	vec4 diff = mat_col * light_col *d; 	\n"
		"	vec4 refl = smoothstep(0.7,0.75,dot(r,light_dir)) * light_col;\n"
		"	return 0.1*refl + 2.25*diff;\n"
		"}\n"
		"\n"
		"void main(void)\n"
		"{\n"
		"	vec4 mat_col = vec4(0.7,0.6,1.0,1.0);\n"
		"	\n"
		"	vec3 light0_dir = vec3(0.0,1.0,0.0);\n"
		"	vec4 light0_col = vec4(0.9,0.95,0.95,1.0);\n"
		"	\n"
		"	vec3 light1_dir = vec3(0.0,0.0,1.0);\n"
		"	vec4 light1_col = vec4(.8,.8,.6,1.0);\n"
		"	\n"
		"	gl_FragColor = \n"
		"	0.5*glazed_shader(mat_col, light0_col, light0_dir)+\n"
		"	0.5*glazed_shader(mat_col, light1_col, light1_dir);\n"
		"	\n"
		"	gl_FragColor.a = 1.0;\n"
		"}\n";
		
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
	glUseProgram(prog);
	
}

void display()
{
	static unsigned int l;
    if(redo_display_list)
    {
        cout << "Creating display list" << endl;
        l = glGenLists(1);
        glNewList(l, GL_COMPILE);
        draw(mesh, per_vertex_normals);
        glEndList();
        redo_display_list = false;
		glutTimerFunc(spin_timer, spin, 0);	
	}
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    ball->set_gl_modelview();
	if(do_wireframe)
	{
		if(GLEW_EXT_geometry_shader4)
			draw_triangles_in_wireframe(mesh, per_vertex_normals, Vec3f(1,0,0));
		else
			draw_wireframe_oldfashioned(mesh, per_vertex_normals, Vec3f(1,0,0));
	}
	else if(!do_textures)
	{
		setupshader();	
		glCallList(l);
	}
	else
		glCallList(l);
	
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
		case '\033': exit(0); break;
		case 'w': do_wireframe = !do_wireframe; break;
		case 'f': per_vertex_normals = !per_vertex_normals; redo_display_list = true; break;
		case 's': 
		{
			ofstream f("ball.out", ios::binary);
			if(f) f.write(reinterpret_cast<const char*>(ball),sizeof(QuatTrackBall));
		}
		break;
		case 'l':
		{
			ifstream f("ball.out", ios::binary);
			if(f) f.read(reinterpret_cast<char*>(ball),sizeof(QuatTrackBall));
		}
			break;
		case 't': do_textures = !do_textures;
			break;
    }
	redo_display_list=true;
}

int main(int argc, char** argv)
{
	Util::ArgExtracter ae(argc, argv);
	
	bool redo_normals = ae.extract("-n");
	
    // GLUT INIT
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
    glutInitWindowSize(win_size_x, win_size_y);
    glutInit(&argc, argv);
    main_window = glutCreateWindow("OBJ Viewer");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(mouse_motion);
    glutMouseFunc(mouse);
    //glutIdleFunc(idle);
	
	glewInit();
	
    // GL INIT
    glClearColor(.8f, 0.9f, 1.0f, 0.f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glShadeModel(GL_SMOOTH);
	
    // LOAD OBJ
    string fn;
    int arg_no = ae.no_remaining_args();
    if(arg_no>1)
        fn = ae.get_last_arg();
    else
        fn = "../../data/head.obj";
	
	load(fn, mesh);
	load_textures(mesh);
	
	if(!mesh.has_normals() || redo_normals)
	{
		cout << "Computing normals" << endl;
		mesh.compute_normals();
	}
	
	// Initialize Trackball
	Vec3f c;
	float r;
	mesh.get_bsphere(c,r);
	r *= 1.5;
	ball = new QuatTrackBall(c,r,800,800);
	
	// Setup projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(53,1.0f,r/100.0,r*3.0);
	glMatrixMode(GL_MODELVIEW);
	
	// Pass control to GLUT
	glutMainLoop();
	
	return 0;
}
