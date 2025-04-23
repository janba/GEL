cmake_minimum_required(VERSION 3.30)

add_executable(covariance_test covariance_test.cpp)
target_include_directories(covariance_test PRIVATE ../..)
set_target_properties(covariance_test PROPERTIES CXX_STANDARD 20)
target_link_libraries(covariance_test GEL)

add_test(NAME GEL.Test.CGLA-covariance COMMAND covariance_test)