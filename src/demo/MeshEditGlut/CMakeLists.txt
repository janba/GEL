find_package (GLEW)

include_directories(../..)
add_executable(glut_main glut_main.cpp)
target_link_libraries(glut_main GEL glut Threads::Threads GLEW::GLEW)
