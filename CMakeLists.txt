cmake_minimum_required(VERSION 3.25) 
project(GEL)

find_package(Threads)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_definitions(-DGLEW_STATIC)
endif ()
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    include_directories(./src/GEL)
    add_definitions(-DGLEW_STATIC)
endif ()
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures for MacOS" FORCE)
endif ()

option(Use_GLGraphics "Compile the OpenGL Viewer" ON)
if (Use_GLGraphics)
    find_package(OpenGL REQUIRED)
    include(FetchContent)
    FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.3.2
    )
    FetchContent_MakeAvailable(glfw)
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        FIND_LIBRARY(COCOA_LIBRARY Cocoa)
        add_compile_definitions(GL_SILENCE_DEPRECATION)
        include_directories(./src/)
        SET(CMAKE_INSTALL_RPATH "@loader_path/")
        # SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    endif ()
endif ()

# Define the GEL library
add_library(GEL STATIC)

# Add include directories for GEL
target_include_directories(GEL PRIVATE ./src)

# Define source lists for GEL library
file(GLOB_RECURSE CGLA_SRC_LIST_CPP ./src/GEL/CGLA/*.cpp)
file(GLOB_RECURSE CGLA_SRC_LIST_C ./src/GEL/CGLA/*.c)
list(APPEND CGLA_SRC_LIST ${CGLA_SRC_LIST_CPP} ${CGLA_SRC_LIST_C})

file(GLOB_RECURSE GEO_SRC_LIST_CPP ./src/GEL/Geometry/*.cpp)
file(GLOB_RECURSE GEO_SRC_LIST_C ./src/GEL/Geometry/*.c)
list(APPEND GEO_SRC_LIST ${GEO_SRC_LIST_CPP} ${GEO_SRC_LIST_C})

file(GLOB_RECURSE MESH_SRC_LIST_CPP ./src/GEL/HMesh/*.cpp)
file(GLOB_RECURSE MESH_SRC_LIST_C ./src/GEL/HMesh/*.c)
list(APPEND MESH_SRC_LIST ${MESH_SRC_LIST_CPP} ${MESH_SRC_LIST_C})

file(GLOB_RECURSE UTIL_SRC_LIST_CPP ./src/GEL/Util/*.cpp)
file(GLOB_RECURSE UTIL_SRC_LIST_C ./src/GEL/Util/*.c)
list(APPEND UTIL_SRC_LIST ${UTIL_SRC_LIST_CPP} ${UTIL_SRC_LIST_C})


# Add source files to GEL library (depending on use of GLGraphics)
if (Use_GLGraphics)
    file(GLOB_RECURSE GL_SRC_LIST_CPP ./src/GEL/GL/*.cpp)
    file(GLOB_RECURSE GL_SRC_LIST_C ./src/GEL/GL/*.c)
    list(APPEND GL_SRC_LIST ${GL_SRC_LIST_CPP} ${GL_SRC_LIST_C})

    file(GLOB_RECURSE GLG_SRC_LIST_CPP ./src/GEL/GLGraphics/*.cpp)
    file(GLOB_RECURSE GLG_SRC_LIST_C ./src/GEL/GLGraphics/*.c)
    list(APPEND GLG_SRC_LIST ${GLG_SRC_LIST_CPP} ${GLG_SRC_LIST_C})
    target_sources(GEL PRIVATE ${GL_SRC_LIST} ${GLG_SRC_LIST} ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MESH_SRC_LIST} ${UTIL_SRC_LIST})
    target_include_directories(GEL PRIVATE ./src/GEL/GL/include)
    target_link_libraries(GEL PRIVATE Threads::Threads OpenGL::GLU OpenGL::GL ${COCOA_LIBRARY})
else ()
    target_sources(GEL PRIVATE ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MESH_SRC_LIST} ${UTIL_SRC_LIST})
    target_link_libraries(GEL PRIVATE Threads::Threads)
endif ()

# Locate the sources for PyGEL
file(GLOB_RECURSE PYG_SRC_LIST ./src/PyGEL/*.cpp)

# Define the PyGEL library. Again two versions depending on use of GLGraphics
if (Use_GLGraphics)
    add_library(PyGEL SHARED ${PYG_SRC_LIST})
    target_include_directories(PyGEL PRIVATE ./src)
    if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        target_include_directories(PyGEL PRIVATE ${glfw_SOURCE_DIR}/include)
    endif ()
    target_link_libraries(PyGEL PRIVATE GEL glfw ${GLFW_LIBRARIES})
else ()
    list(REMOVE_ITEM PYG_SRC_LIST "./src/PyGEL/Viewer.cpp")
    add_library(PyGEL SHARED ${PYG_SRC_LIST})
    target_include_directories(PyGEL PRIVATE ./src)
    target_link_libraries(PyGEL PRIVATE GEL)
endif ()

# Install GEL library
install(TARGETS GEL LIBRARY DESTINATION lib)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/GEL
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    )

# Add test projects
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(src/test)
endif()


# No install instructions for PyGEL. It is left in the build directory and
# included in the wheel package.
