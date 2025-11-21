cmake_minimum_required(VERSION 3.30)

get_filename_component(GEL_INCLUDE_DIR ../.. ABSOLUTE)

add_executable(RsR_test RsR_test.cpp)

set_target_properties(RsR_test PROPERTIES
    CXX_STANDARD 20
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined"
)
target_include_directories(RsR_test PRIVATE "${GEL_INCLUDE_DIR}" "${DOCTEST_INCLUDE_DIR}")

target_link_libraries(RsR_test PRIVATE GEL nanobench)

doctest_discover_tests(RsR_test TEST_PREFIX "RsR.")
#add_test(NAME GEL.Test.RsR COMMAND RsR_test)