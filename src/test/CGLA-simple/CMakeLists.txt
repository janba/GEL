cmake_minimum_required(VERSION 3.30)
project(CGLA-covarianceGEL)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)

# Locate the GEL library
find_library(GEL_LIB GEL PATHS ../../../../build)

add_executable(simple_test simple_test.cpp)
target_link_libraries(simple_test ${GEL_LIB})
