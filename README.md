## INTRODUCTION
GEL is a C++ library of geometry processing tools intended for computer graphics applications. In particular, GEL
contains a fairly mature half-edge based library, an efficient kD tree and data structures for volumetric data.
Functionality includes quadric based simplification, mesh optimization, distance field computation, iso surface polygonization 
and more. A linear algebra library for small vectors and matrices is also included as well as tools for visualizing meshes
using OpenGL.

PyGEL3D is a set of Python bindings for a subset of the features in GEL. In particular, PyGEL covers almost all the mesh features. In addition PyGEL has its own viewer based on OpenGL and PyGEL can be used from Jupyter notebooks. In this case, it is possible to visualize meshes using a plotly widget. A significant benefit here is that when the notebook is exported to HTML, the 3D view comes along.

## DOCUMENTATION
Some installation instructions below. But for more documentation please see the doc directory. There is a doxygen script for creating a reference manual and a latex file intro.tex which explains the basics. Please doxygen or pdflatex your documentation. A license is also found in the intro document.

### Building on XCode/OSX
An XCode project has been created and is found in the GEL_MAC directory. The XCode project produces a framework for the GEL library. There are also Xcode projects for the demos and for PyGEL. However, if you want to install PyGEL as a wheel package, you shoud use cmake as discussed below.

### Building on Windows
An Visual Studio 2013 project has been created and is found in the GEL_WIN directory. The project produces a framework for the GEL library.
Look in the doc folder for installation of glut and glew

### Building on Linux
The CMake files (see below) work cross platform.

### Building with CMake
There is a basic CMakeLists.txt - it compiles the GEL library and the PyGEL library but not the demos.
```
#$> mkdir build; cd build; cmake ..; make -j8 ; cd ..
```
### Creating a PyGEL3D package and installing it
to create a wheel package for PyGEL3D that works with Python 3, you can next issue the command
```
#$> python3 setup.py bdist_wheel
```
Install using something like
```
#$> sudo -H pip install dist/PyGEL3D-*.whl
```
For this to work, you need to have wheel and setuptools installed. 
## Requirements
Compiling both GEL and PyGEL requires that you have OpenGL and GLFW installed unless you choose not to compile graphics support which you can do by setting `Use_GLGraphics` to `OFF`.

PyGEL3D has a module called 'js' which produces graphics suitable for Jupyter notebooks. This module is based on plotly which must then be installed for it to work. You will also need numpy. However, these required libraries will be downloaded when you install PyGEL3D if you use pip.
