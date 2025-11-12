# Spatial Module

::: pygel3d.spatial

The `spatial` module provides spatial data structures for efficient geometric queries.

## I3DTree Class

The `I3DTree` class is a kD-tree specialized for mapping 3D points to integers (typically array indices).

### Creating a Tree

```python
from pygel3d import I3DTree
import numpy as np

# Create tree
tree = I3DTree()

# Insert points with associated indices
points = np.random.rand(1000, 3)
for i, point in enumerate(points):
    tree.insert(point, i)
```

## Tree Operations

### Insertion
- `tree.insert(point, value)` - Insert a point with an associated integer value

### Queries
- `tree.closest_point(query_point)` - Find nearest point, returns associated value
- `tree.in_sphere(center, radius)` - Find all points within radius

### Information
- `tree.size()` - Number of points in tree

## MeshDistance Class

The `MeshDistance` class enables efficient distance computations to triangle meshes.

### Creating a Distance Object

```python
from pygel3d import MeshDistance
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Create distance object
dist = MeshDistance(m)
```

## Distance Queries

### Point Queries
- `dist.signed_distance(point)` - Signed distance to mesh (negative inside)
- `dist.distance(point)` - Unsigned distance to mesh
- `dist.ray_inside_test(point, direction)` - Test if ray from point hits mesh

### Information
- Returns scalar distance values
- Signed distance uses mesh orientation (closed meshes only)

## Example Usage

### k-Nearest Neighbors

```python
from pygel3d import I3DTree
import numpy as np

# Generate random points
points = np.random.rand(1000, 3) * 10

# Build tree
tree = I3DTree()
for i, p in enumerate(points):
    tree.insert(p, i)

# Find nearest point to query
query = [5.0, 5.0, 5.0]
nearest_idx = tree.closest_point(query)
nearest_point = points[nearest_idx]

print(f"Nearest point to {query}: {nearest_point}")
print(f"Index: {nearest_idx}")
```

### Points in Sphere

```python
from pygel3d import I3DTree
import numpy as np

# Build tree
tree = I3DTree()
points = np.random.rand(1000, 3) * 10
for i, p in enumerate(points):
    tree.insert(p, i)

# Find all points within radius
center = [5.0, 5.0, 5.0]
radius = 2.0
indices = tree.in_sphere(center, radius)

print(f"Found {len(indices)} points within radius {radius}")
```

### Distance Field Computation

```python
from pygel3d import MeshDistance
import pygel3d.hmesh as hmesh
import numpy as np

# Load mesh
m = hmesh.load("bunny.obj")

# Create distance object
dist = MeshDistance(m)

# Create grid
grid_size = 50
x = np.linspace(-1, 1, grid_size)
y = np.linspace(-1, 1, grid_size)
z = np.linspace(-1, 1, grid_size)

# Compute distance field
distances = np.zeros((grid_size, grid_size, grid_size))
for i, xi in enumerate(x):
    for j, yj in enumerate(y):
        for k, zk in enumerate(z):
            point = [xi, yj, zk]
            distances[i, j, k] = dist.signed_distance(point)

print(f"Distance field shape: {distances.shape}")
print(f"Min distance: {distances.min()}")
print(f"Max distance: {distances.max()}")
```

### Inside/Outside Testing

```python
from pygel3d import MeshDistance
import pygel3d.hmesh as hmesh

# Load closed mesh
m = hmesh.load("sphere.obj")
dist = MeshDistance(m)

# Test points
test_points = [
    [0, 0, 0],      # Inside
    [10, 10, 10],   # Outside
    [0.5, 0.5, 0.5] # Near surface
]

for point in test_points:
    d = dist.signed_distance(point)
    status = "inside" if d < 0 else "outside"
    print(f"Point {point}: {status} (distance: {abs(d):.4f})")
```

## Performance Tips

- **Tree Construction**: Build the tree once and reuse for multiple queries
- **Batch Queries**: Use vectorized operations when possible
- **Mesh Simplification**: For distance queries, simpler meshes are faster
- **Precision**: Consider using lower precision for large-scale computations
