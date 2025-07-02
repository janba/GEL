cmake_minimum_required(VERSION 3.30)

add_executable(Assert_test Assert_test.cpp)
target_include_directories(Assert_test PRIVATE ../.. ${DOCTEST_INCLUDE_DIR})
set_target_properties(Assert_test PROPERTIES
        CXX_STANDARD 20)
target_link_libraries(Assert_test GEL nanobench)

doctest_discover_tests(Assert_test TEST_PREFIX "Assert.")