cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(collapse_test collapse_test.cpp)

set_target_properties(collapse_test PROPERTIES
    CXX_STANDARD 20
)

target_include_directories(collapse_test PRIVATE ${GEL_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR})
target_link_libraries(collapse_test GEL)

doctest_discover_tests(collapse_test TEST_PREFIX "Hierarchical.Collapse.")
