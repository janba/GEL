cmake_minimum_required(VERSION 3.30)
project(OGL-testGEL)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)
find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)

add_executable(ogl_test ogl_test.cpp)
target_link_libraries(ogl_test ${GEL_LIB} OpenGL::GL OpenGL::GLU GLUT::GLUT)
