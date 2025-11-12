# Quick Start Guide

This guide will help you get started with PyGEL3D through practical examples.

## Basic Mesh Operations

### Loading and Saving Meshes

```python
import pygel3d.hmesh as hmesh

# Load a mesh from file
mesh = hmesh.obj_load("bunny.obj")

# You can also use the generic load function
# which detects the file format automatically
mesh = hmesh.load("model.obj")

# Save the mesh
hmesh.obj_save("output.obj", mesh)
```

Supported formats include OBJ, OFF, PLY, and X3D.

### Creating a Simple Mesh

```python
import pygel3d.hmesh as hmesh

# Create a new empty mesh
m = hmesh.Manifold()

# Add vertices (returns vertex IDs)
v0 = m.add_vertex([0, 0, 0])
v1 = m.add_vertex([1, 0, 0])
v2 = m.add_vertex([0.5, 1, 0])

# Add a face using vertex positions
face_id = m.add_face([0, 0, 0,  # v0
                      1, 0, 0,  # v1
                      0.5, 1, 0])  # v2

print(f"Created mesh with {m.no_vertices()} vertices and {m.no_faces()} faces")
```

### Mesh Information

```python
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Get basic statistics
print(f"Vertices: {m.no_vertices()}")
print(f"Faces: {m.no_faces()}")
print(f"Halfedges: {m.no_halfedges()}")

# Check mesh properties
print(f"Valid: {hmesh.valid(m)}")
print(f"Closed: {hmesh.closed(m)}")
print(f"Boundary curves: {hmesh.count_boundary_curves(m)}")

# Get bounding box
bbox_min, bbox_max = hmesh.bbox(m)
print(f"Bounding box: {bbox_min} to {bbox_max}")
```

## Mesh Processing

### Smoothing

```python
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Catmull-Clark smoothing
hmesh.cc_smooth(m)

# Laplacian smoothing
hmesh.laplacian_smooth(m, weight=0.5, iter=10)

# Taubin smoothing (better volume preservation)
hmesh.taubin_smooth(m, iter=10)
```

### Subdivision

```python
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Catmull-Clark subdivision
hmesh.cc_split(m)

# Loop subdivision (for triangle meshes)
hmesh.loop_split(m)

# Root-3 subdivision
hmesh.root3_subdivide(m)
```

### Simplification

```python
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Quadric error metric simplification
# Keep 50% of the original faces
hmesh.quadric_simplify(m, keep_fraction=0.5)
```

### Triangulation

```python
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Triangulate the mesh
hmesh.shortest_edge_triangulate(m)

# Alternative triangulation method
hmesh.ear_clip_triangulate(m)
```

## Visualization

### OpenGL Viewer

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl

# Load a mesh
m = hmesh.load("bunny.obj")

# Create a viewer and display
viewer = gl.Viewer()
viewer.display(m, mode='g')  # 'g' for glazed (shaded)
```

Display modes:
- `'w'`: Wireframe
- `'n'`: Normal (flat shading)
- `'g'`: Glazed (smooth shading)
- `'i'`: Isophote lines
- `'l'`: Line field
- `'s'`: Scalar field

### Jupyter Notebook Visualization

```python
import pygel3d.hmesh as hmesh
import pygel3d.jupyter_display as jd

# Load a mesh
m = hmesh.load("bunny.obj")

# Display in Jupyter
jd.display(m, smooth=True)
```

## Graph Processing

### Creating and Working with Graphs

```python
import pygel3d.graph as graph
import pygel3d.hmesh as hmesh

# Create a graph from a mesh skeleton
m = hmesh.load("model.obj")
g = graph.from_mesh(m)

# Access graph properties
print(f"Nodes: {g.no_nodes()}")
print(f"Edges: {g.no_edges()}")

# Get node positions
positions = g.positions()

# Save and load graphs
graph.save("skeleton.graph", g)
g2 = graph.load("skeleton.graph")
```

### Graph to Mesh Conversion

```python
import pygel3d.graph as graph
import pygel3d.hmesh as hmesh

# Load a graph
g = graph.load("skeleton.graph")

# Convert to a cylindrical mesh
m = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, m, fudge=0.5)

# Save the result
hmesh.save("output.obj", m)
```

## Spatial Queries

### Distance Computation

```python
import pygel3d.hmesh as hmesh
from pygel3d import MeshDistance
import numpy as np

# Load a mesh
m = hmesh.load("model.obj")

# Create a distance object
dist = MeshDistance(m)

# Query distance from a point
point = [0, 0, 0]
distance = dist.signed_distance(point)

print(f"Distance from origin: {distance}")
```

### kD-Tree Queries

```python
from pygel3d import I3DTree
import numpy as np

# Create some 3D points
points = np.random.rand(1000, 3)

# Build a kD-tree
tree = I3DTree()
for i, p in enumerate(points):
    tree.insert(p, i)

# Query nearest point
query_point = [0.5, 0.5, 0.5]
nearest_idx = tree.closest_point(query_point)

print(f"Nearest point index: {nearest_idx}")
```

## Next Steps

- Explore the [API Reference](../api/overview.md) for detailed function documentation
- Check out the [Tutorials](../tutorials/mesh-operations.md) for in-depth guides
- See [Examples](../examples.md) for complete working examples
