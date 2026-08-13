cmake_minimum_required(VERSION 3.25)
project(HRsR)
set(CMAKE_CXX_STANDARD 20)
include_directories(../..)

add_executable(HRsR HRsR.cpp)

# Prefer linking against an in-tree GEL target when this demo is added
# from the main project; otherwise link explicitly to the local archive.
if(TARGET GEL)
	target_link_libraries(HRsR PRIVATE GEL)
else()
	set(GEL_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/libGEL.a")
	if(NOT EXISTS "${GEL_LIB_PATH}")
		message(FATAL_ERROR "Could not find GEL static library at ${GEL_LIB_PATH}. Build GEL first with cmake --build build.")
	endif()
	target_link_libraries(HRsR PRIVATE "${GEL_LIB_PATH}")
endif()
