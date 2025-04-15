cmake_minimum_required(VERSION 3.30)
project(GELtest)

# Locate the GEL library if invoked from here
if(NOT GEL)
    find_library(GEL PATHS ../../../build GEL_PATH)
endif()

add_subdirectory(CGLA-covariance)
#add_subdirectory(CGLA-mat)
#add_subdirectory(CGLA-ogl)
#add_subdirectory(CGLA-simple)
#add_subdirectory(CGLA-vec)
#add_subdirectory(Geometry-kdtree)
#add_subdirectory(GLGraphics-console)
