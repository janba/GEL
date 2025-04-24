cmake_minimum_required(VERSION 3.30)

add_executable(mat_test mat_test.cpp)
target_include_directories(mat_test PRIVATE ../..)
set_target_properties(covariance_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(mat_test GEL)

add_test(NAME GEL.Test.CGLA-mat COMMAND mat_test)
