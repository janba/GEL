# PyGEL3D Documentation

Welcome to the PyGEL3D documentation. PyGEL3D is a Python library for geometry processing, providing tools for working with 3D polygonal meshes, graphs, and spatial data structures.

## What is PyGEL3D?

PyGEL3D is a collection of classes and functions for geometry processing tasks, especially those involving 3D polygonal meshes. It is based on the C++ GEL library and provides a Python interface for most of its functionality.

## Key Features

GEL/PyGEL3D contains several data structures for spatial data. Our mature halfedge data structure for polygonal meshes and a graph data structure are the most important components. 

Several algorithms that build on these data structure are offered by the library. Some of these are based directly on our own research. For instance, the Rotations System Reconstruction (RSR) algorithm for reconstruction of triangle meshes from point clouds and the Local Separator Skeletonization method for extracting curve skeletons from anything represented as a spatially embedded graph are both implemented in this library.

A number of standard techniques have also been implemented. For instance, Garland-Heckbert simplification, several schemes for subdivision, curvature computation methods, and a number of methods for smoothing (including feature preserving smoothing) are also provided by the library.

Finally, both in the C++ and Python library, tools are provided for analysis and editing a mesh. Since the half-edge based representation is not restricted to triangles, these tools work for general polygonal meshes.

## Modules

PyGEL3D consists of five main modules:

- **[hmesh](api/hmesh.md)**: Halfedge mesh data structure and operations
- **[graph](api/graph.md)**: Graph data structure for spatial graphs
- **[spatial](api/spatial.md)**: Spatial data structures (kD-Trees, distance queries)
- **[gl_display](api/gl_display.md)**: OpenGL-based visualization
- **[jupyter_display](api/jupyter_display.md)**: Jupyter notebook integration

## Quick Example

```python
from pygel3d import hmesh, gl_display as gl
from sys import argv

m = hmesh.load(argv[1])

hmesh.close_holes(m)
hmesh.triangulate(m)
hmesh.quadric_simplify(m, 0.05)

v = gl.Viewer()
v.display(m) # Hit ESC to exit
```

## Getting Started

Ready to dive in? Check out the [Installation Guide](getting-started/installation.md) and [Quick Start Tutorial](getting-started/quickstart.md).

## Links

- [GitHub Repository](https://github.com/janba/GEL)
- [PyPI Package](https://pypi.org/project/PyGEL3D/)
