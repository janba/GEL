project(GEL)
cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 14)

find_package (OpenGL)

if (APPLE)
find_package (glfw3 REQUIRED)
	FIND_LIBRARY(COCOA_LIBRARY Cocoa)
	add_compile_definitions(GL_SILENCE_DEPRECATION)
	include_directories(./src/)
endif (APPLE)

if (UNIX)
find_package (glfw3)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif (UNIX)

if(WIN32)
	include_directories(./src/GEL)
	add_compile_definitions(-DGLEW_STATIC)
endif()

include_directories(./src/GEL/GL/include)

aux_source_directory(./src/GEL/GL GL_SRC_LIST)
aux_source_directory(./src/GEL/CGLA CGLA_SRC_LIST)
aux_source_directory(./src/GEL/Geometry GEO_SRC_LIST)
aux_source_directory(./src/GEL/HMesh MSH_SRC_LIST)
aux_source_directory(./src/GEL/GLGraphics GLG_SRC_LIST)
aux_source_directory(./src/GEL/Util UTL_SRC_LIST)

if(WIN32)
	add_library(GEL STATIC ${GL_SRC_LIST} ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${GLG_SRC_LIST} ${UTL_SRC_LIST})
	target_link_libraries(GEL OpenGL::GLU OpenGL::GL ${COCOA_LIBRARY})
else()
	add_library(GEL SHARED ${GL_SRC_LIST} ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${GLG_SRC_LIST} ${UTL_SRC_LIST})
	target_link_libraries(GEL OpenGL::GLU OpenGL::GL ${COCOA_LIBRARY})
endif()

# not right now ...
# add_subdirectory(./src/demo)

include_directories(./src)
aux_source_directory(./src/PyGEL PYG_SRC_LIST)

if(WIN32)
	add_library(PyGEL SHARED ${PYG_SRC_LIST})
	target_link_libraries(PyGEL GEL libglfw3)
else()
	add_library(PyGEL SHARED ${PYG_SRC_LIST})
	target_link_libraries(PyGEL GEL glfw)
endif()
