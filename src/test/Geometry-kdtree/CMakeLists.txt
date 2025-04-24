cmake_minimum_required(VERSION 3.30)

add_executable(kdtree-test kdtree-test.cpp)
target_include_directories(kdtree-test PRIVATE ../..)
set_target_properties(kdtree-test PROPERTIES CXX_STANDARD 20)
target_link_libraries(kdtree-test GEL)

add_test(NAME GEL.Test.Geometry-kdtree COMMAND kdtree-test)
