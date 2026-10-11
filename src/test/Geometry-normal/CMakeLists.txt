cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(normal_test normal_test.cpp)

set_target_properties(normal_test PROPERTIES
    CXX_STANDARD 20
)

target_include_directories(normal_test PRIVATE ${GEL_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR})
target_link_libraries(normal_test GEL)

doctest_discover_tests(normal_test TEST_PREFIX "Geometry.Normal.")
