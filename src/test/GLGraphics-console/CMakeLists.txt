cmake_minimum_required(VERSION 3.30)
project(GLGraphics-consoleGEL)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)
find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)

add_executable(console_test console_test.cpp)
target_link_libraries(console_test ${GEL_LIB} OpenGL::GL OpenGL::GLU GLUT::GLUT)
