cmake_minimum_required(VERSION 3.30)

add_executable(vec_test vec_test.cpp)
target_include_directories(vec_test PRIVATE ../..)
set_target_properties(vec_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(vec_test GEL)

add_test(NAME GEL.Test.CGLA-vec COMMAND vec_test)

