cmake_minimum_required(VERSION 3.30)
project(MeshEditGlutGEL)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)
find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)
find_package (Threads REQUIRED)

# Add all .cpp files in the directory to the executable
file(GLOB SOURCES "*.cpp")
add_executable(pathtracer ${SOURCES})
target_link_libraries(pathtracer ${GEL_LIB} OpenGL::GL OpenGL::GLU GLUT::GLUT Threads::Threads)
