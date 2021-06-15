## INTRODUCTION
GEL is a C++ library of geometry processing tools mostly intended for computer graphics applications. In particular, GEL contains a fairly mature half-edge based polygonal mesh (sub)library, a graph data structure, and various spatial data structures.

Functionality includes a number of mesh algorithms such as Garland-Heckbert simplification, mesh optimization, distance field computation, iso-surface polygonization, a function for computing curve skeletons, and more. 

A linear algebra library for small vectors and matrices (2D, 3D, and 4D) is also included as well as tools for visualizing meshes using OpenGL. 

GEL requires a recent C++ compiler but has very few dependencies. For visualization, OpenGL and GLFW are required, but these dependencies can be omitted if visualization is not needed. This somewhat limits the capabilities of GEL, and most projects that use GEL would also require other libraries (such as a linear algebra library) but by not requiring these to be installed as dependencies of GEL, compilation is simplified and porting to new platforms becomes much easier.
### PyGEL
PyGEL is a set of Python bindings for a _subset_ of the features in GEL. In particular, PyGEL covers almost all the mesh features. In addition PyGEL has its own viewer based on OpenGL, and PyGEL can be used from Jupyter notebooks. In this case, it is possible to visualize meshes using a plotly widget. A significant benefit here is that when the notebook is exported to HTML, the 3D view comes along. In fact, this makes PyGEL a useful tool for teaching geometry processing since work can be carried out in a Jupyter notebook and assignments submitted as HTML files. Moreover, PyGEL works in Google [colab](https://colab.research.google.com) notebooks.

## DOCUMENTATION
Some installation instructions below. But for more documentation please see the doc directory. There is a doxygen script for creating a reference manual and a latex file intro.tex which explains the basics. Please doxygen or pdflatex your documentation. A license is also found in the intro document.

## INSTALLING PYGEL WITH PIP

PyGEL is on PyPI and can be installed with pip. For most potential users, there is no reason to look further than:
```
pip install PyGEL3D
```
The PyPI package is called PyGEL3D and not PyGEL. Turns out just PyGEL was taken by a completely unrelated library. 
## BUILDING AND INSTALLING GEL AND PYGEL
If you need or wish to compile GEL there are a few options:
- For all platforms, there are the CMake files. These can be used to build both GEL and PyGEL but not the demos or tests. The setup.py script which is used to make a PyPI package is only useful with this option for building PyGEL.
- For Mac OS there is an Xcode project (in GEL_MAC). There are also projects for demos, tests, and PyGEL.
- For Windows there is a Visual Studio Solution (in GEL_WIN) for GEL and demos.
- For Unices there is a python script (in GEL_UNIX) for building GEL and PyGEL

Especially if you intend to be a user of GEL or PyGEL and not a developer, I _strongly_ recommend that you use CMake to build. GEL and PyGEL are built every commit using GitHub's continuous integration based on the CMake files. Thus, this should be a robust way to compile the project.

Since GEL is primarily developed on Mac OS, the Xcode projects are actively maintained and cover everything. **The two last build options are not actively maintained**. However, they may prove of some use and hence not removed.

### Building with CMake
If you are using a unix-like command line, build with
```
mkdir build; cd build; cmake ..; make -j 8 ; sudo make install; cd ..
```
### Creating a PyGEL3D package and installing it
You can next issue the command
```
python setup.py bdist_wheel
```
Install using something like
```
pip install dist/PyGEL3D-*.whl
```
For this to work, you need to have wheel and setuptools installed. Also `python` and `pip` should be version 3 variants.
## Practical Issues
Compiling both GEL and PyGEL requires that you have OpenGL installed unless you choose not to compile graphics support which you can do by setting `Use_GLGraphics` to `OFF` in the CMake file. GLFW is also needed, but CMake fetches GLFW from github and compiles it along with the GEL code. If you compile in some of the other ways (e.g. using XCode, Visual Studio) there is no simple way to avoid the dependency on graphics libraries. Thus, if you need to avoid the OpenGL requirements, CMake is the way to go.

GEL comes with a few demo applications. In addition to the requirements above, several of these also require GLUT to be installed. Going forward, we should remove the GLUT dependency and move to GLFW for the applications.

PyGEL has a module called `jupyter_display` which produces graphics suitable for Jupyter notebooks. This module is based on plotly which must then be installed for it to work. You will also need numpy. However, if you use pip, these required libraries will be downloaded automatically when you install PyGEL.