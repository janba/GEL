cmake_minimum_required(VERSION 3.30)


add_executable(OBJParser_test OBJParser_test.cpp)
target_include_directories(OBJParser_test PRIVATE ../.. ${DOCTEST_INCLUDE_DIR})
set_target_properties(OBJParser_test PROPERTIES
    CXX_STANDARD 20)
target_link_libraries(OBJParser_test GEL nanobench)

doctest_discover_tests(OBJParser_test TEST_PREFIX "ObjParser.")