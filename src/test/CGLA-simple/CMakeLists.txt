cmake_minimum_required(VERSION 3.30)

add_executable(simple_test simple_test.cpp)
target_include_directories(simple_test PRIVATE ../..)
set_target_properties(simple_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(simple_test GEL)

add_test(NAME GEL.Test.CGLA-simple COMMAND simple_test)
