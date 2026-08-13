cmake_minimum_required(VERSION 3.25)
project(GEL_Demos LANGUAGES CXX)

# Convenience project that builds every C++ demo against an installed GEL
# package. Each subdirectory is also a standalone CMake project:
#
#   cmake -S src/demo/MeshDistance -B src/demo/MeshDistance/build
#
# GEL installs to ~/.local by default. HINTS makes find_package succeed
# without CMAKE_PREFIX_PATH after a standard install.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GEL REQUIRED HINTS "$ENV{HOME}/.local")

add_subdirectory(HRsR)
add_subdirectory(MeshDistance)
add_subdirectory(MeshEditGlut)
add_subdirectory(OBJViewer)
add_subdirectory(Pathtracer)
add_subdirectory(Raytracer)
add_subdirectory(RotationSystemReconstruction)
add_subdirectory(Skeletonization)
