cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(cycle_test cycle_test.cpp)

set_target_properties(cycle_test PROPERTIES
    CXX_STANDARD 20
)

target_include_directories(cycle_test PRIVATE ${GEL_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR})
target_link_libraries(cycle_test GEL)

doctest_discover_tests(cycle_test TEST_PREFIX "Geometry.Cycles.")
