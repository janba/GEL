# Mesh Operations Tutorial

This tutorial covers common mesh operations in detail.

## Loading Meshes

PyGEL3D supports several mesh formats:

```python
import pygel3d.hmesh as hmesh

# Auto-detect format
m = hmesh.load("model.obj")

# Specific format loaders
m = hmesh.obj_load("model.obj")
m = hmesh.off_load("model.off")
m = hmesh.ply_load("model.ply")
m = hmesh.x3d_load("model.x3d")
```

## Mesh Traversal

### Iterating Over Elements

```python
# Iterate over vertices
for v in m.vertices():
    pos = m.positions()[3*v:3*v+3]
    print(f"Vertex {v}: {pos}")

# Iterate over faces
for f in m.faces():
    print(f"Face {f}")

# Iterate over halfedges
for h in m.halfedges():
    print(f"Halfedge {h}")
```

### Circulating Around Elements

```python
# Circulate around a vertex
vertex_id = 0
neighbors = m.circulate_vertex(vertex_id, mode='v')
print(f"Neighbor vertices: {neighbors}")

# Circulate around a face
face_id = 0
vertices = m.circulate_face(face_id, mode='v')
print(f"Face vertices: {vertices}")
```

## Mesh Repair

### Stitching Coincident Vertices

```python
# Merge vertices that are very close
threshold = 1e-6
num_merged = hmesh.stitch_mesh(m, threshold)
print(f"Merged {num_merged} vertex pairs")
```

### Closing Holes

```python
# Close holes up to 100 edges
max_hole_size = 100
hmesh.close_holes(m, max_hole_size)
```

### Removing Degenerate Features

```python
# Remove caps (nearly flat triangles at edges)
hmesh.remove_caps(m, threshold=0.1)

# Remove needles (very thin triangles)
hmesh.remove_needles(m, threshold=0.1, averagePositions=True)
```

## Mesh Smoothing

### Catmull-Clark Smoothing

```python
# Apply CC smoothing (subdivides and smooths)
hmesh.cc_smooth(m)
```

### Laplacian Smoothing

```python
# Smooth with weight 0.5 for 10 iterations
hmesh.laplacian_smooth(m, weight=0.5, iter=10)
```

### Taubin Smoothing

```python
# Better volume preservation
hmesh.taubin_smooth(m, iter=10)
```

## Mesh Subdivision

### Catmull-Clark

```python
# Subdivide using Catmull-Clark scheme
hmesh.cc_split(m)
```

### Loop Subdivision

```python
# For triangle meshes
hmesh.loop_split(m)
```

## Mesh Simplification

```python
# Simplify to 50% of faces
hmesh.quadric_simplify(
    m,
    keep_fraction=0.5,
    singular_thresh=1e-10,
    error_thresh=1e10
)
```

## Mesh Optimization

### Edge Flipping

```python
# Optimize valency (tries to make vertices degree 6)
hmesh.optimize_valency(m, anneal=True)

# Maximize minimum angle
hmesh.maximize_min_angle(m, thresh=0.1, anneal=True)

# Minimize curvature
hmesh.minimize_curvature(m, anneal=True)
```

See the [API Reference](../api/hmesh.md) for complete function documentation.
