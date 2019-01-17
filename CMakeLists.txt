project(GEL)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 14)

find_package (OpenGL)
find_package (glfw3)


if (APPLE)
FIND_LIBRARY(COCOA_LIBRARY Cocoa)
add_compile_definitions(GL_SILENCE_DEPRECATION)
include_directories(./src/)
endif (APPLE)

if (UNIX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif (UNIX)

if(WIN32)
	include_directories(./src/GEL)
    # TODO for Windows
endif()

include_directories(./src/GEL/GL/include)

aux_source_directory(./src/GEL/GL GL_SRC_LIST)
aux_source_directory(./src/GEL/CGLA CGLA_SRC_LIST)
aux_source_directory(./src/GEL/Geometry GEO_SRC_LIST)
aux_source_directory(./src/GEL/HMesh MSH_SRC_LIST)
aux_source_directory(./src/GEL/GLGraphics GLG_SRC_LIST)
aux_source_directory(./src/GEL/Util UTL_SRC_LIST)

add_library(GEL SHARED ${GL_SRC_LIST} ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${GLG_SRC_LIST} ${UTL_SRC_LIST})
target_link_libraries(GEL OpenGL::GLU OpenGL::GL ${COCOA_LIBRARY})

#add_subdirectory(./src/demo)

include_directories(./src)
aux_source_directory(./src/PyGEL PYG_SRC_LIST)
add_library(PyGEL SHARED ${PYG_SRC_LIST})
target_link_libraries(PyGEL GEL glfw)
