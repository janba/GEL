find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)
find_package (GLEW REQUIRED)

include_directories(../..)
add_executable(glut_main glut_main.cpp)
target_link_libraries(glut_main GEL OpenGL::GL OpenGL::GLU GLUT::GLUT Threads::Threads GLEW::GLEW)
