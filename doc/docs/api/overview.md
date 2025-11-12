# API Overview

PyGEL3D provides five main modules for geometry processing:

## Module Summary

### [hmesh](hmesh.md) - Halfedge Mesh
The core mesh data structure and operations. This module provides the `Manifold` class for representing polygonal meshes using the halfedge data structure, along with numerous functions for mesh manipulation, analysis, and processing.

**Key Classes:**
- `Manifold`: Halfedge-based polygonal mesh
- `MeshDistance`: Compute distances to triangle meshes

**Key Functions:**
- Mesh I/O (load, save in various formats)
- Mesh processing (smoothing, subdivision, simplification)
- Mesh queries (topology, geometry, measurements)
- Mesh modification (editing, optimization)

### [graph](graph.md) - Graph Processing
Spatial graph data structure and algorithms. Useful for representing curve skeletons and other graph-based geometric structures.

**Key Classes:**
- `Graph`: 3D spatial graph with nodes and edges

**Key Functions:**
- Graph construction from meshes
- Skeletonization algorithms
- Graph-to-mesh conversion
- Graph optimization

### [spatial](spatial.md) - Spatial Data Structures
Efficient spatial queries and data structures.

**Key Classes:**
- `I3DTree`: kD-tree for 3D point-to-integer mapping
- `MeshDistance`: Distance queries to triangle meshes

**Key Functions:**
- Nearest neighbor queries
- Distance computations
- Spatial indexing

### [gl_display](gl_display.md) - OpenGL Visualization
Interactive 3D visualization using OpenGL.

**Key Classes:**
- `Viewer`: OpenGL-based mesh and graph viewer

**Key Functions:**
- Interactive mesh display
- Multiple rendering modes
- Camera control

### [jupyter_display](jupyter_display.md) - Jupyter Integration
Visualization tools for Jupyter notebooks using Plotly.

**Key Functions:**
- Interactive 3D widgets in notebooks
- Export notebooks with embedded 3D graphics
- Compatible with Google Colab

## Common Workflows

### Mesh Processing Pipeline
```python
import pygel3d.hmesh as hmesh

# 1. Load mesh
m = hmesh.load("input.obj")

# 2. Clean up
hmesh.stitch_mesh(m, 1e-6)
hmesh.close_holes(m)

# 3. Process
hmesh.cc_smooth(m)
hmesh.triangulate(m)

# 4. Optimize
hmesh.quadric_simplify(m, keep_fraction=0.5)

# 5. Save
hmesh.save("output.obj", m)
```

### Graph-Based Skeletonization
```python
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph

# 1. Load mesh
m = hmesh.load("model.obj")

# 2. Extract skeleton
g = graph.from_mesh(m)

# 3. Process graph
graph.smooth(g, iter=10)
graph.prune(g)

# 4. Convert back to mesh
result = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, result)
```

### Distance Field Computation
```python
import pygel3d.hmesh as hmesh
from pygel3d import MeshDistance
import numpy as np

# Load mesh
m = hmesh.load("model.obj")

# Create distance object
dist = MeshDistance(m)

# Query distances
points = np.random.rand(1000, 3) * 10
distances = [dist.signed_distance(p) for p in points]
```

## Data Types

### Common Types
- **Vertex/Face/Halfedge IDs**: Represented as integers (size_t)
- **3D Positions**: Lists or tuples of 3 floats [x, y, z]
- **Vectors**: Lists of values
- **Bounding Boxes**: Tuple of two 3D positions (min, max)

### Index-Based Access
Most PyGEL3D functions use index-based access:
```python
# Get vertex position by index
vertex_id = 0
pos = m.positions()[vertex_id * 3:(vertex_id + 1) * 3]

# Iterate over all vertices
for v in m.vertices():
    # v is a vertex index
    pass
```

## Performance Considerations

- **In-place Operations**: Most mesh operations modify the mesh in-place
- **Memory**: Large meshes may require significant memory
- **Python Overhead**: For performance-critical code, consider working with the C++ library directly
- **Vectorization**: When possible, use batch operations rather than loops

## Error Handling

PyGEL3D functions generally:
- Return `False` or `-1` on failure for operations that can fail
- Raise exceptions for invalid inputs
- Print warnings for non-critical issues

Always check return values:
```python
success = hmesh.obj_load("file.obj", m)
if not success:
    print("Failed to load mesh")
```

## Next Steps

Explore each module's detailed documentation:
- [HMesh Module](hmesh.md)
- [Graph Module](graph.md)
- [Spatial Module](spatial.md)
- [GL Display Module](gl_display.md)
- [Jupyter Display Module](jupyter_display.md)
