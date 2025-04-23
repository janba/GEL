cmake_minimum_required(VERSION 3.30)

find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

add_executable(console_test console_test.cpp)
target_include_directories(console_test PRIVATE ../..)
set_target_properties(console_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(console_test GEL OpenGL::GL OpenGL::GLU GLUT::GLUT)

add_test(NAME GEL.Test.GLGraphics-console COMMAND console_test)
