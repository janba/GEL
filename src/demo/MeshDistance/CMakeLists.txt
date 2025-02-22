cmake_minimum_required(VERSION 3.30)
project(MeshDistanceGEL)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)

add_executable(meshdist meshdist.cpp)
target_link_libraries(meshdist ${GEL_LIB})
