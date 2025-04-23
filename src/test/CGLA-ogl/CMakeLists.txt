cmake_minimum_required(VERSION 3.30)

find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

add_executable(ogl_test ogl_test.cpp)
target_include_directories(ogl_test PRIVATE ../..)
set_target_properties(ogl_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(ogl_test GEL OpenGL::GL OpenGL::GLU GLUT::GLUT)

add_test(NAME GEL.Test.CGLA-OGL COMMAND ogl_test)