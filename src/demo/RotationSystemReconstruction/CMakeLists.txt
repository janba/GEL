cmake_minimum_required(VERSION 3.25)
project(RotationSystemReconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GEL REQUIRED HINTS "$ENV{HOME}/.local")

add_executable(rsr rsr.cpp)
target_link_libraries(rsr PRIVATE GEL::GEL)
