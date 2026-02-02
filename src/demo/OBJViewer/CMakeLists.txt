cmake_minimum_required(VERSION 3.30)
project(OBJViewerGEL)
set(CMAKE_FIND_FRAMEWORK NEVER)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)
find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)
find_package (Threads REQUIRED)

add_executable(objviewer objviewer.cpp)
target_link_libraries(objviewer ${GEL_LIB} OpenGL::GL OpenGL::GLU GLUT::GLUT Threads::Threads)
