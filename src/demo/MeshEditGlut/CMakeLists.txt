find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)
find_package (GLEW REQUIRED)
find_package (Threads REQUIRED)


include_directories(../..)
link_directories(../../../../build)

add_executable(glut_main glut_main.cpp)
target_link_libraries(glut_main GEL OpenGL::GL OpenGL::GLU GLUT::GLUT Threads::Threads GLEW::GLEW)

install(TARGETS glut_main RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
