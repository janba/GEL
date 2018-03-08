project(GEL)
cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_STANDARD 14)

find_package (Threads)
find_package (OpenGL)
find_package (GLEW REQUIRED)

if (APPLE)
FIND_LIBRARY(COCOA_LIBRARY Cocoa)
endif (APPLE)

if (UNIX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lpthread -Wall")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -Wall")
endif (UNIX)
if(WIN32)
    # TODO for Windows
endif()

include_directories(./src/GEL)
include_directories(./src/GEL/GL/include)

aux_source_directory(./src/GEL/CGLA CGLA_SRC_LIST)
aux_source_directory(./src/GEL/Geometry GEO_SRC_LIST)
aux_source_directory(./src/GEL/HMesh MSH_SRC_LIST)
aux_source_directory(./src/GEL/GLGraphics GLG_SRC_LIST)
aux_source_directory(./src/GEL/Util UTL_SRC_LIST)

add_library(GEL SHARED ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${GLG_SRC_LIST} ${UTL_SRC_LIST})# ./src/GEL/GL/src/glad.c)
target_link_libraries(GEL Threads::Threads OpenGL::GLU OpenGL::GL GLEW::GLEW ${COCOA_LIBRARY})

add_subdirectory(./src/demo)

include_directories(./src)
aux_source_directory(./src/PyGEL PYG_SRC_LIST)
add_library(PyGEL SHARED ${PYG_SRC_LIST})
target_link_libraries(PyGEL GEL)

