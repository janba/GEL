cmake_minimum_required(VERSION 3.25)
project(HRsR LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GEL REQUIRED HINTS "$ENV{HOME}/.local")

add_executable(HRsR HRsR.cpp)
target_link_libraries(HRsR PRIVATE GEL::GEL)
