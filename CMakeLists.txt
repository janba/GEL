cmake_minimum_required(VERSION 3.11)
project(GEL)
find_package(Threads)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_CXX_FLAGS "-std=c++1z")
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    include_directories(./src/GEL)
    add_definitions(-DGLEW_STATIC)
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")
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
    FetchContent_GetProperties(glfw)
    if(NOT glfw_POPULATED)
        FetchContent_Populate(glfw)

        set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
        set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
        set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
        set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

        add_subdirectory(${glfw_SOURCE_DIR} ${glfw_BINARY_DIR})
    endif()
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        FIND_LIBRARY(COCOA_LIBRARY Cocoa)
        add_compile_definitions(GL_SILENCE_DEPRECATION)
        include_directories(./src/)
        SET(CMAKE_INSTALL_RPATH "@loader_path/")
        SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    endif ()
endif ()

aux_source_directory(./src/GEL/CGLA CGLA_SRC_LIST)
aux_source_directory(./src/GEL/Geometry GEO_SRC_LIST)
aux_source_directory(./src/GEL/HMesh MSH_SRC_LIST)
aux_source_directory(./src/GEL/Util UTL_SRC_LIST)

if (Use_GLGraphics)
    include_directories(./src/GEL/GL/include)
    aux_source_directory(./src/GEL/GL GL_SRC_LIST)
    aux_source_directory(./src/GEL/GLGraphics GLG_SRC_LIST)
    add_library(GEL STATIC ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${UTL_SRC_LIST} ${GL_SRC_LIST} ${GLG_SRC_LIST})
    target_link_libraries(GEL Threads::Threads OpenGL::GLU OpenGL::GL ${COCOA_LIBRARY})
    # not right now ...
    # add_subdirectory(./src/demo)
else ()
    add_library(GEL STATIC ${CGLA_SRC_LIST} ${GEO_SRC_LIST} ${MSH_SRC_LIST} ${UTL_SRC_LIST})
    target_link_libraries(GEL Threads::Threads)
endif ()


target_include_directories(GEL PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)
aux_source_directory(./src/PyGEL PYG_SRC_LIST)

if (Use_GLGraphics)
    add_library(PyGEL SHARED ${PYG_SRC_LIST})
    if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        target_include_directories(PyGEL PRIVATE ${glfw_SOURCE_DIR}/include)
    endif ()
    target_link_libraries(PyGEL GEL glfw ${GLFW_LIBRARIES})
else ()
    list(REMOVE_ITEM PYG_SRC_LIST "./src/PyGEL/Viewer.cpp")
    add_library(PyGEL SHARED ${PYG_SRC_LIST})
    target_link_libraries(PyGEL GEL)
endif ()

target_include_directories(PyGEL PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)
    
install(TARGETS GEL)

install(TARGETS PyGEL GEL
        LIBRARY DESTINATION PyGEL
        ARCHIVE DESTINATION PyGEL)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/GEL
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    )
