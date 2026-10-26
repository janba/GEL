cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(solve_test solve_test.cpp)

set_target_properties(solve_test PROPERTIES
    CXX_STANDARD 20
)

target_include_directories(solve_test PRIVATE ${GEL_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR})
target_link_libraries(solve_test GEL)

doctest_discover_tests(solve_test TEST_PREFIX "CGLA.Solve.")
