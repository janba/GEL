cmake_minimum_required(VERSION 3.25)
if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING
        "Minimum macOS version" FORCE)
endif()
project(GEL VERSION 0.7.3 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Default to a user-local prefix so headers and libraries land in
# ~/.local/{include,lib} without sudo. Override with -DCMAKE_INSTALL_PREFIX.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT WIN32
        AND DEFINED ENV{HOME} AND NOT "$ENV{HOME}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH
        "Install path prefix, prepended onto install directories." FORCE)
endif()

find_package(Threads REQUIRED)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64" CACHE STRING "Build architectures for MacOS" FORCE)
endif ()

option(Use_GLGraphics "Compile the OpenGL Viewer" ON)
if (Use_GLGraphics)
    find_package(OpenGL REQUIRED)
    include(FetchContent)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.3.10
    )
    FetchContent_MakeAvailable(glfw)
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        SET(CMAKE_INSTALL_RPATH "@loader_path/")
        # SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    endif ()
endif ()

# Define the GEL library
add_library(GEL STATIC)
add_library(GEL::GEL ALIAS GEL)

target_compile_features(GEL PUBLIC cxx_std_20)

# Consumers include headers as <GEL/HMesh/Manifold.h>, so the include root
# is src/ in the build tree and ${CMAKE_INSTALL_INCLUDEDIR} after install.
target_include_directories(GEL
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

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
    target_link_libraries(GEL PUBLIC Threads::Threads OpenGL::GLU OpenGL::GL)
    if (APPLE)
        target_link_libraries(GEL PUBLIC "-framework Cocoa")
        target_compile_definitions(GEL PUBLIC GL_SILENCE_DEPRECATION)
    endif ()
    if (NOT APPLE)
        # Match the previous Linux/Windows GLEW build; consumers that include
        # <GEL/GL/glew.h> need the same definition.
        target_compile_definitions(GEL PUBLIC GLEW_STATIC)
    endif ()
else ()
    target_sources(GEL PRIVATE ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MESH_SRC_LIST} ${UTIL_SRC_LIST})
    target_link_libraries(GEL PUBLIC Threads::Threads)
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
    list(FILTER PYG_SRC_LIST EXCLUDE REGEX ".*/Viewer\\.cpp$")
    add_library(PyGEL SHARED ${PYG_SRC_LIST})
    target_include_directories(PyGEL PRIVATE ./src)
    target_link_libraries(PyGEL PRIVATE GEL)
endif ()

# TODO: this one needs to be specifically modified so we do not install it globally, but rather
# TODO: inside our own directories. It should probably live in include/GEL/Util/phmap
## Include the parallel-hashmap library
#include(FetchContent)
## Our headers depend on it, so we need to it alongside GEL
#set(PHMAP_INSTALL TRUE CACHE INTERNAL "Install parallel_hashmap")
#FetchContent_Declare(
#    parallel-hashmap
#    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
#    GIT_TAG        v2.0.0    # adjust tag/branch/commit as needed
#)
#FetchContent_MakeAvailable(parallel-hashmap)
## Pass the directory to GEL during a build
#target_include_directories(
#    GEL PUBLIC
#    $<BUILD_INTERFACE:${parallel-hashmap_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
#    )


# ---------------------------------------------------------------------------
# Install GEL as a CMake package (find_package(GEL))
# ---------------------------------------------------------------------------
set(GEL_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/GEL")

install(TARGETS GEL
    EXPORT GELTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/GEL
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
    )

install(EXPORT GELTargets
    FILE GELTargets.cmake
    NAMESPACE GEL::
    DESTINATION ${GEL_INSTALL_CMAKEDIR}
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GELConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GELConfig.cmake"
    INSTALL_DESTINATION ${GEL_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GELConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GELConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GELConfigVersion.cmake"
    DESTINATION ${GEL_INSTALL_CMAKEDIR}
)

# Also export a build-tree package so find_package(GEL) works with
# -DCMAKE_PREFIX_PATH=<build-dir> before a full install.
export(EXPORT GELTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/GELTargets.cmake"
    NAMESPACE GEL::
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
