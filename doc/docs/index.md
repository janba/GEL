# PyGEL3D Documentation

Welcome to the PyGEL3D documentation. PyGEL3D is a Python library for geometry processing, providing tools for working with 3D polygonal meshes, graphs, and spatial data structures.

## What is PyGEL3D?

PyGEL3D is a collection of classes and functions for geometry processing tasks, especially those involving 3D polygonal meshes. It is based on the C++ GEL library and provides a Python interface for most of its functionality.

## Key Features

### **Halfedge Mesh Data Structure**
- Rich API for mesh traversal and manipulation
- Support for polygonal meshes (not just triangles)
- Extensive mesh processing operations

### **Graph Processing**
- 3D spatial graphs with vertices and edges
- Local Separator Skeletonization for curve skeleton extraction
- Inverse skeletonization (Face Extrusion Quad meshes from graphs)

### **Spatial Data Structures**
- kD-Tree for efficient spatial queries
- Distance computations to triangle meshes
- Bounding volume hierarchies

### **Visualization**
- OpenGL-based viewer for interactive 3D visualization
- Jupyter notebook integration with Plotly widgets
- Export notebooks to HTML with embedded 3D graphics

### **Advanced Features**
- Garland-Heckbert mesh simplification
- Signed distance field computation
- Iso-surface polygonization
- Mesh smoothing (anisotropic, Taubin, and more)
- Topological analysis

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
