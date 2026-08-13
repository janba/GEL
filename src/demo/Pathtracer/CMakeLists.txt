cmake_minimum_required(VERSION 3.25)
project(Pathtracer LANGUAGES CXX)

if(APPLE)
	set(CMAKE_FIND_FRAMEWORK FIRST)
else()
	set(CMAKE_FIND_FRAMEWORK NEVER)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GEL REQUIRED HINTS "$ENV{HOME}/.local")
find_package(OpenGL REQUIRED)

if(APPLE)
	find_library(APPLE_GLUT_FRAMEWORK GLUT REQUIRED)
	find_library(APPLE_COREFOUNDATION_FRAMEWORK CoreFoundation REQUIRED)
else()
	find_package(GLUT REQUIRED)
endif()

file(GLOB SOURCES "*.cpp")
add_executable(pathtracer ${SOURCES})

if(APPLE)
	target_link_libraries(pathtracer PRIVATE
		GEL::GEL
		OpenGL::GL
		OpenGL::GLU
		${APPLE_GLUT_FRAMEWORK}
		${APPLE_COREFOUNDATION_FRAMEWORK})
else()
	target_link_libraries(pathtracer PRIVATE
		GEL::GEL
		OpenGL::GL
		OpenGL::GLU
		GLUT::GLUT)
endif()
