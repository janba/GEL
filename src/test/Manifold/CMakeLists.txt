cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(Manifold_test Manifold_test.cpp)

set_target_properties(Manifold_test PROPERTIES
    CXX_STANDARD 20
)

target_include_directories(Manifold_test PRIVATE ${GEL_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR})
target_link_libraries(Manifold_test GEL nanobench)

doctest_discover_tests(Manifold_test TEST_PREFIX "Manifold.")
