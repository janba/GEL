cmake_minimum_required(VERSION 3.30)

add_executable(ThreadPool_test ThreadPool_test.cpp)
target_include_directories(ThreadPool_test PRIVATE ../.. ${DOCTEST_INCLUDE_DIR})
set_target_properties(ThreadPool_test PROPERTIES
        CXX_STANDARD 20)
target_link_libraries(ThreadPool_test GEL nanobench)

#add_test(NAME GEL.Test.ThreadPool COMMAND ThreadPool_test)
doctest_discover_tests(ThreadPool_test TEST_PREFIX "ThreadPool.")