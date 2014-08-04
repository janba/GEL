#include <typeinfo>
#include <iostream>

#include <GEL/CGLA/Mat4x4f.h>
#include <GEL/CGLA/Mat2x2f.h>

#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec2i.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3f.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace CGLA;


Mat4x4f perspective_Mat4x4f(float d)
{
  Mat4x4f m(0.0f);
  
  /* Eye at the origin, looking down the negative z axis */

  m[0][0] = 1.0;
  m[1][1] = 1.0;
  m[2][2] = 1.0;
  m[3][2] = -1.0/d;
   
  return m;
}


static void display( void )
{

  // ----------------------------------------
  // 0. Set up viewing parameters

  Vec3f up(0,1,0);             // The direction that is most nearly up ..
  Vec3f eye(3,3,3);            // position of eye 
  Vec3f centre(0,0,0);         // what we are looking at 
  float image_plane_dist = 1;  // distance from eye to image plane.


  // ----------------------------------------
  // 1. Create view coordinate system
  //
  // Note that the args in the cross(.,.) call
  // do not commute

  Vec3f n = centre - eye;
  n.normalize();

  Vec3f u = cross(n,up);
  u.normalize();

  Vec3f v = cross(u,n);

  //----------------------------------------
  // 2. Create matrices

  // Create viewing matrix. We use the basis change method.
  // Notice how the direction of z is flipped. That is because
  // we look down the -z direction
  Mat4x4f mview(Vec4f(u,0), Vec4f(v,0), Vec4f(-n,0), Vec4f(0,0,0,1));

  //Create translation matrix. 
  Mat4x4f mtrans = translation_Mat4x4f(centre-eye);

  // Create projection matrix
  Mat4x4f mpers  = perspective_Mat4x4f(image_plane_dist);

  // Concatenate the translation, viewing and projection matrices
  Mat4x4f m = mpers * mview * mtrans;

	std::cout << mview << mtrans << mpers << m << std::endl;

  //----------------------------------------
  // 3. Create points 

	Vec4f axes[3] = { Vec4f(2, 0, 0,1), Vec4f(0, 2, 0, 1), Vec4f(0, 0, 2, 1) };
  Vec4f paxes[3];
  Vec4f p[9] = { Vec4f(0, 0, 0, 1), Vec4f(1, 0, 0,1), Vec4f(0, 1, 0, 1),
		Vec4f(1,1,0,1), Vec4f(0,0,1,1), Vec4f(1,0,1,1),  
		Vec4f(0,1,1,1), Vec4f(1,1,1,1)};
  Vec4f pp[9];

  //----------------------------------------
  // 4. project and dehomogenize points
  
  paxes[0] = m * axes[0];
  paxes[1] = m * axes[1];
  paxes[2] = m * axes[2];
  paxes[0].de_homogenize();
  paxes[1].de_homogenize();
  paxes[2].de_homogenize();

  for (int i=0;i<9;i++) 
    {
	  pp[i] = m * p[i];
      pp[i].de_homogenize();
    }


  //----------------------------------------
  // 5. Draw _projected_ points in 2D using OpenGL

  // Clear screen
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glBegin(GL_LINES);
  glColor3f(1,0,0);
  glVertex2fv(pp[0].get());
  glVertex2fv(paxes[0].get());

  glColor3f(0,1,0);
  glVertex2fv(pp[0].get());
  glVertex2fv(paxes[1].get());
  
  glColor3f(0,0,1);
  glVertex2fv(pp[0].get());
  glVertex2fv(paxes[2].get());

  glColor3f(0,0,0);  
  for(int i=0;i<4;i++)
    {
      glVertex2fv(pp[2*i           + 0 ].get());
      glVertex2fv(pp[2*i           + 1 ].get());
    }
  for(int i=0;i<4;i++)
    {
      glVertex2fv(pp[(4*(i/2) + i%2) + 0 ].get());
      glVertex2fv(pp[(4*(i/2) + i%2) + 2 ].get());
    }
  for(int i=0;i<4;i++)
    { 
      glVertex2fv(pp[1*i           + 0 ].get());
      glVertex2fv(pp[1*i           + 4 ].get());
    }
  glEnd();
	glFlush();
}

static void reshape( int width, int height )
{
  glViewport( 0, 0, width, height );
}


static void key( unsigned char key, int x, int y )
{
  switch (key) {
  case 27:
    exit(0);
    break;
  }
}
static void init_GL()
{
  // Initialize GL, i.e. setup projection
  // and possibly other things as well
  glClearColor(1,1,1,1);
  gluOrtho2D(-3,3,-3,3);    
}
static void init_GLUT(int argc, char *argv[])
{
  // Initialize glut, open and create window.
  glutInit( &argc, argv );
  glutInitWindowSize( 400 , 400 );
  glutCreateWindow(argv[0]);

  // Register callback functions
  glutReshapeFunc( reshape );
  glutKeyboardFunc( key );
  glutDisplayFunc( display );
}

int main( int argc, char *argv[] )
{
  init_GLUT(argc, argv);
  init_GL();
  glutMainLoop();
  return 0;
}





