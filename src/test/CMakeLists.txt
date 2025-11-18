cmake_minimum_required(VERSION 3.30)

include(FetchContent)

# Add benchmarking library
FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG v4.1.0
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

# Add testing library
FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.12
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

FetchContent_MakeAvailable(nanobench)
# # Add nanobench to a project with:
# target_link_libraries(some_executable nanobench)

FetchContent_MakeAvailable(doctest)
FetchContent_GetProperties(doctest SOURCE_DIR doctest_source_dir)
set(DOCTEST_INCLUDE_DIR ${doctest_source_dir}/doctest CACHE INTERNAL "Path to include folder for doctest")
# # Add doctest to a project with:
# target_include_directories(some_executable PRIVATE ${DOCTEST_INCLUDE_DIR})
include(${doctest_source_dir}/scripts/cmake/doctest.cmake)

add_subdirectory(CGLA-covariance)
add_subdirectory(CGLA-mat)
add_subdirectory(CGLA-simple)
add_subdirectory(CGLA-vec)
add_subdirectory(Geometry-kdtree)
add_subdirectory(Util-Assert)
add_subdirectory(Manifold)
# The following two rely on key input in order to exit
#add_subdirectory(CGLA-ogl)
#add_subdirectory(GLGraphics-console)

add_subdirectory(Util-ThreadPool)
add_subdirectory(Util)
