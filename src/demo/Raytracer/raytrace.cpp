#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <GEL/GLGraphics/GLViewController.h>
#include <GEL/GL/glew.h>
#ifdef __APPLE__
#include <GLUT/GLUT.h>
#else
#include <GL/glut.h>
#endif

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat4x4f.h>

#include <GEL/HMesh/Manifold.h>

#include <GEL/Geometry/TriMesh.h>
#include <GEL/Geometry/obj_load.h>
#include <GEL/Geometry/Ray.h>
#include <GEL/Geometry/BSPTree.h>
#include <GEL/Geometry/build_bbtree.h>
#include <GEL/Geometry/AABox.h>
#include <GEL/Geometry/OBox.h>

#include <GEL/Util/Timer.h>

#include "Camera.h"

//#define USE_BSP

using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace HMesh;
using namespace GLGraphics;

/*
 * TODO:
 * - BBOX remove HMesh dependency - that is crazy.
 * - BBox visit child nodes in order of how far away the intersection point on
 *        the bbox is. Closest nodes visited first. Cull nodes farther than
 *        an actual intersection point.
 * - BBox Smooth interpolation of triangle normals. Straightforward.
 */


namespace
{
    const int MAX_OBJECTS = 4;   // Maximum number of triangles in a BSP tree node
    const int MAX_LEVEL = 20;    // Maximum number of BSP tree subdivisions
    
    const unsigned int TEX_SIZE = 512;
    
    bool raytrace = false;
    bool done = false;
    bool shadow = false;
    
    unsigned int winx = TEX_SIZE;     // Screen width
    unsigned int winy = TEX_SIZE;     // Screen height
    
    unsigned int PIXEL_SUBDIVS = 1;
    
    int mouse_state = GLUT_UP;
    int mouse_button = 0;
    int spin_timer = 20;
    
    GLViewController *vctrl;
    
    TriMesh mesh;
    vector<const TriMesh*> mesh_vector(1, &mesh);
    vector<Mat4x4f> transforms(1, identity_Mat4x4f());
    
    double light_pow = 1.0;
    Vec3f light_dir = normalize(Vec3f(1.0, 1.0, 1.0));
    Vec3d background(0.8, 0.9, 1.0);
    
    BSPTree tree;
    Camera* cam;
    
    Vec3f image[TEX_SIZE][TEX_SIZE];
    unsigned int image_tex;
    
    // Function to generate a random number between 0 and 1.
    // gel_rand() returns a random number ranging from 0 to GEL_RAND_MAX.
    inline double my_random()
    {
        return gel_rand()/static_cast<double>(GEL_RAND_MAX);
    }
    
    AABBTree bb_tree;
}

void spin(int x);

//////////////////////////////////////////////////////////////
//      I N I T I A L I Z A T I O N
//////////////////////////////////////////////////////////////

void init_texture(unsigned int& tex)
{
    glGenTextures(1, reinterpret_cast<GLuint*>(&tex));
    
    glBindTexture(GL_TEXTURE_2D, tex);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    
    // load the texture image
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
	             TEX_SIZE, TEX_SIZE,
	             0, GL_RGB, GL_FLOAT, image[0][0].get());
}

void initRaytracer()
{
    Vec3f c;
    float r;
    mesh.get_bsphere(c, r);
    r *= 1.5;
    
    // Initialize track ball
    vctrl = new GLViewController(winx, winy, c, r);
    
    // Initialize corresponding camera
    cam = new Camera(c - Vec3f(r), c, Vec3f(0.0f, 1.0f, 0.0f), 1.0f);
    
#ifdef USE_BSP
    cout << "Constructing BSP tree..." << endl;
    tree.init(mesh_vector, transforms, MAX_OBJECTS, MAX_LEVEL);
    tree.build();
#else
	// AABB TREE
	Manifold m;
	vector<int> faces(mesh.geometry.no_faces(), 3);
	cout << "Creating manifold" << endl;
    m.build(mesh.geometry.no_vertices(),
            reinterpret_cast<const float*>(&mesh.geometry.vertex(0)),
            faces.size(),
            &faces[0],
            reinterpret_cast<const int*>(&mesh.geometry.face(0)));
	cout << "Building tree" << endl;
	build_AABBTree(m, bb_tree);
#endif
}

void initGL()
{
    glShadeModel(GL_SMOOTH);
    glDisable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glColor3f(0.0, 0.0, 0.0);
}


//////////////////////////////////////////////////////////////
//      S H A D E   F U N C T I O N S
//////////////////////////////////////////////////////////////

double shadow_shade(Ray& r)
{
    r.compute_position();
    Ray shadow(r.hit_pos, light_dir);
    double s = tree.intersect(shadow) ? 0.0 : 1.0;
    
    r.compute_normal();
    return s*light_pow*dot(r.hit_normal, light_dir);
}

double lambertian_shade(Ray& r)
{
#ifdef USE_BSP
	r.compute_normal();
#endif
    return light_pow*dot(r.hit_normal, light_dir);
}

double (*shade_ray[2])(Ray&) = { lambertian_shade,
    shadow_shade      };

//////////////////////////////////////////////////////////////
//      D R A W   F U N C T I O N S
//////////////////////////////////////////////////////////////

/*
 void enable_textures(TriMesh& tm)
 {
 for(unsigned int i=0;i<tm.materials.size(); ++i)
 {
 Material& mat = tm.materials[i];
 if(mat.tex_name != "")
 {
 string name = mat.tex_path + mat.tex_name;
 
 GLuint tex_id;
 if(load_image_into_texture(name, tex_id))
 mat.tex_id = tex_id;
 }
 }
 }
 */

void set_perspective_proj()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    gluPerspective(64.0, 1.0, 0.1, 1000.0);
    
    glMatrixMode(GL_MODELVIEW);
}

void set_ortho_proj()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    
    glMatrixMode(GL_MODELVIEW);
}

void draw_texture(unsigned int tex)
{
    static GLfloat verts[] = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
    
    glColor4f(1.0, 1.0, 1.0, 1.0);
    
    glBindTexture(GL_TEXTURE_2D, tex);
    glEnable(GL_TEXTURE_2D);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    
    glVertexPointer(2, GL_FLOAT, 0, verts);
    glTexCoordPointer(2, GL_FLOAT, 0, verts);
    
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    
    glDisable(GL_TEXTURE_2D);
}

void drawOBJ()
{
    static bool washere = false;
    static unsigned int disp_list;
    
    if(!washere)
    {
        disp_list = glGenLists(1);
        glNewList(disp_list, GL_COMPILE);
        
        glBegin(GL_TRIANGLES);
        for(int i = 0; i < mesh.geometry.no_faces(); ++i)
        {
            Vec3i n_face = mesh.normals.face(i);
            Vec3i g_face = mesh.geometry.face(i);
            for(int j=0;j<3;j++)
            {
                double shade = 0.5;
                
                if(n_face != Geometry::NULL_FACE)
                {
                    Vec3f norm = normalize(mesh.normals.vertex(n_face[j]));
                    glNormal3fv(norm.get());
                    shade = light_pow*dot(norm, light_dir);
                }
                
                glColor3d(shade, shade, shade);
                Vec3f vert = mesh.geometry.vertex(g_face[j]);
                glVertex3fv(vert.get());
            }
        }
        glEnd();
        
        glEndList();
        washere = true;
    }
    glCallList(disp_list);
}


//////////////////////////////////////////////////////////////
//      G L U T   C A L L B A C K   F U N C T I O N S
//////////////////////////////////////////////////////////////

void display()
{
    static bool first = true;
    
    if(first)
    {
        first = false;
        glutTimerFunc(spin_timer, spin, 0);
    }
    
    Vec3f eye, focus, up;
	vctrl->get_view_param(eye, focus, up);
	cam->set(eye, focus, up, cam->get_focal_dist());
	
    if(raytrace)
    {
        raytrace = false;
        
        cout << "Raytracing";
        float win_to_vp = 1.0f/static_cast<float>(TEX_SIZE);
        float lowerleft = 0.5 + win_to_vp*0.5;
        float step = win_to_vp/static_cast<float>(PIXEL_SUBDIVS);
        
        vector<Vec2f> jitter(PIXEL_SUBDIVS*PIXEL_SUBDIVS);
        for(unsigned int i = 0; i < PIXEL_SUBDIVS; ++i)
            for(unsigned int j = 0; j < PIXEL_SUBDIVS; ++j)
            {
                jitter[i*PIXEL_SUBDIVS + j][0] = (my_random() + (j%PIXEL_SUBDIVS))*step;
                jitter[i*PIXEL_SUBDIVS + j][1] = (my_random() + (i%PIXEL_SUBDIVS))*step;
            }
        
        Util::Timer tim;
        tim.start();
        for(unsigned int i = 0; i < TEX_SIZE; ++i)
        {
            for(unsigned int j = 0; j < TEX_SIZE; ++j)
            {
                Vec3d sum(0.0f);
                Vec2f vp_pos(j*win_to_vp - lowerleft, i*win_to_vp - lowerleft);
                
                for(unsigned int ky = 0; ky < PIXEL_SUBDIVS; ++ky)
                    for(unsigned int kx = 0; kx < PIXEL_SUBDIVS; ++kx)
                    {
                        Ray r = cam->get_ray(vp_pos + jitter[ky*PIXEL_SUBDIVS + kx]);
                        
#ifdef USE_BSP
                        if(tree.intersect(r))
                            sum += Vec3d(shade_ray[shadow](r));
                        else
                            sum += background;
#else
                        float t = FLT_MAX;
                        bb_tree.intersect(r);
                        if(r.has_hit)
                            sum += Vec3d(shade_ray[0](r));
                        else
                            sum += background;
#endif
                    }
                image[i][j] = Vec3f(sum/static_cast<double>(PIXEL_SUBDIVS*PIXEL_SUBDIVS));
            }
            if(((i + 1) % 50) == 0) cerr << ".";
        }
        cout << " - " << tim.get_secs() << " secs " << endl;
        cout << endl;
        
        init_texture(image_tex);
        
        done = true;
    }
    
    if(done)
    {
        set_ortho_proj();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        
        draw_texture(image_tex);
    }
    else
    {
        glEnable(GL_DEPTH_TEST);
        
        cam->glSetPerspective(winx, winy);
        
        glClearColor(background[0], background[1], background[2], 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        
        cam->glSetCamera();
        
        glColor3f(0.5, 0.5, 0.5);
        drawOBJ();
        
        glDisable(GL_DEPTH_TEST);
    }
    
    glutSwapBuffers();
}

void reshape(int w, int h)
{
    winx = w; winy = h;
    
    vctrl->reshape(winx, winy);
    
    glViewport(0, 0, winx, winy);
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
        case '+':
            ++PIXEL_SUBDIVS;
            cout << "Rays per pixel: " << PIXEL_SUBDIVS*PIXEL_SUBDIVS << endl;
            break;
        case '-':
            if(PIXEL_SUBDIVS > 1)
                --PIXEL_SUBDIVS;
            cout << "Rays per pixel: " << PIXEL_SUBDIVS*PIXEL_SUBDIVS << endl;
            break;
        case 'r':
            if(done) done = false;
            else raytrace = true;
            break;
        case 's':
            shadow = !shadow;
            cout << "Shadow " << (shadow ? "on." : "off.") << endl;
            break;
        case 27:
            delete vctrl;
            delete cam;
            exit(0);
    }
}

void mouse(int btn, int state, int x, int y)
{
    if(state == GLUT_DOWN)
    {
        if(btn == GLUT_LEFT_BUTTON)
            vctrl->grab_ball(ROTATE_ACTION, Vec2i(x, y));
        else if(btn == GLUT_MIDDLE_BUTTON)
            vctrl->grab_ball(ZOOM_ACTION, Vec2i(x, y));
        else if(btn == GLUT_RIGHT_BUTTON)
            vctrl->grab_ball(PAN_ACTION, Vec2i(x, y));
    }
    else if(state == GLUT_UP)
        vctrl->release_ball();
    
    mouse_state = state;
    mouse_button = btn;
    
    glutPostRedisplay();
}

void move(int x, int y)
{
    vctrl->roll_ball(Vec2i(x, y));
    
    glutPostRedisplay();
}

void spin(int x)
{
    vctrl->try_spin();
    glutTimerFunc(spin_timer, spin, 0);
    glutPostRedisplay();
}


//////////////////////////////////////////////////////////////
//                        M A I N
//////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
#ifdef __BORLANDC__
    _control87(MCW_EM, MCW_EM);  // Borland C++ will crash OpenGL if this
    // magic line is not inserted
#endif
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(winx, winy);
    glutCreateWindow("Press 'r' to raytrace");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(move);
    
    // LOAD OBJ
    string filename;
    if(argc > 1)
    {
        filename = argv[1];
        cout << "Loading " << filename << endl;
        
        obj_load(filename, mesh);
        //if(!mesh.has_normals())
        //{
        cout << "Computing normals" << endl;
        mesh.compute_normals();
        //}
        
        cout << "No. of triangles: " << mesh.geometry.no_faces() << endl;
    }
    else
    {
        obj_load("../../data/dolphins.obj", mesh);
        
        cout << "Computing normals" << endl;
        mesh.compute_normals();
        
        //cout << "Usage: raytrace any_object.obj";
        //exit(0);
    }
    
    initRaytracer();
    initGL();    
    
    glutMainLoop();
    
    return 0;
}

