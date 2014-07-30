#include <GEL/GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <GEL/GLGraphics/Console.h>
#include <GEL/CGLA/Vec3f.h>

GLGraphics::Console console;

int width, height;
bool console_visible = true;

static void display()
{
    glClearColor(0.4f, 0.4f, 0.4f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (console_visible)
    {
        //draw console in upper half of screen
        glPushAttrib(GL_VIEWPORT_BIT);
        glViewport(0, height-height/2,
                   width, height/2);
        console.display();
        glPopAttrib();
    }

    assert(glGetError() == GL_NO_ERROR);
    glutSwapBuffers();
}

static void reshape(int w, int h)
{
    width = w;
    height = h;

    glViewport( 0, 0, width, height);
}

static void keyboard(unsigned char key, int x, int y)
{
    //toggle console with ESC
    if (key == 27)
    {
        console_visible = !console_visible;
        glutPostRedisplay();
        return;
    }

    if (console_visible)
    {
        console.keyboard(key);
        glutPostRedisplay();
        return;
    }

    //switch (key)
    //{
    //}

    glutPostRedisplay();
}

static void special(int key, int x, int y)
{
    if (console_visible)
        switch (key) {
            case GLUT_KEY_UP:
                console.key_up();
                break;
            case GLUT_KEY_DOWN:
                console.key_down();
                break;
            case GLUT_KEY_LEFT:
                console.key_left();
                break;
            case GLUT_KEY_RIGHT:
                console.key_right();
                break;
            case GLUT_KEY_HOME:
                console.key_home();
                break;
            case GLUT_KEY_END:
                console.key_end();
                break;
                
            default:
                break;
        }
    glutPostRedisplay();
}

void vararg_test(const std::vector<std::string>& args)
{
    console.printf("Number of arguments: %i", int(args.size()));
    for (int i=0; i<int(args.size()); ++i)
        console.printf("  arg %i: %s", i, args[i].c_str());
}

int main( int argc, char *argv[] )
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(768, 768);
  glutInitWindowPosition(256, 256);
  glutCreateWindow(argv[0]);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special);
  glutDisplayFunc(display);

  console.printf("GLGraphics console test.");
  console.newline();

  //some examples of use:

  console.reg_cmd0("quit",
      std::bind(&std::exit, EXIT_SUCCESS), "Exit application.");

  console.reg_cmd1<int>("quit",
      std::bind(&std::exit, std::placeholders::_1),
      "Exit application with specified exit code.");

  console.reg_cmd0("fullscreen", std::bind(&glutFullScreen), "Switch to fullscreen.");

#if 1
  //needs lambda
  console.reg_cmd0("window_pos", [&]{
      console.printf("window_pos = %i %i", glutGet(GLUT_WINDOW_X),
          glutGet(GLUT_WINDOW_Y));
  }, "Show window position.");
  console.reg_cmd0("window_size", [&]{
      console.printf("window_size = %i %i", glutGet(GLUT_WINDOW_WIDTH),
          glutGet(GLUT_WINDOW_HEIGHT));
  }, "Show window position.");
#endif

  console.reg_cmd2<int,int>("window_pos", std::bind(&glutPositionWindow,
      std::placeholders::_1, std::placeholders::_2),
      "Set the window position.");

  console.executef("window_pos %i %i", 384, 256);

  console.reg_cmd2<int,int>("window_size", std::bind(&glutReshapeWindow,
      std::placeholders::_1, std::placeholders::_2),
      "Set the window size.");

  console.reg_cmdN("vararg_test", vararg_test, "Test of variable number of arguments.");

  using namespace GLGraphics;

  Console::variable<int> test_int(42);
  test_int.reg(console, "test_int", "Some clever help string..");
  console.execute("test_int");
  console.execute("test_int 167");

  Console::variable<float> test_float(3.14f);
  test_float.reg(console, "test_float", "Well..");
  console.execute("test_float");
  console.execute("test_float 2.71");

  Console::variable<std::string> test_string("Hello, world!");
  test_string.reg(console, "test_string", "Well..");
  console.execute("test_string");
  console.execute("test_string \"some other string with spaces in\"");

  Console::variable<CGLA::Vec3f> test_Vec3f(CGLA::Vec3f(1,2,3));
  test_Vec3f.reg(console, "test_Vec3f", "Well..");
  console.execute("test_Vec3f");
  console.execute("test_Vec3f [ 12 1234 15]");

  CGLA::Vec3f in(0,1,2);
  std::stringstream ss;
  ss << in;
  CGLA::Vec3f out;
  ss >> out;
  assert(in == out);

  glutMainLoop();
  return 0;
}





