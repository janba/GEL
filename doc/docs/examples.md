# Examples

This page provides complete, runnable examples demonstrating common PyGEL3D workflows.

## Basic Mesh Processing

### Load, Process, and Save

```python
import pygel3d.hmesh as hmesh
from sys import argv

def print_stats(m):
    # Print statistics
    print(f"Vertices: {len(m.vertices())}")
    print(f"Faces: {len(m.faces())}")
    print(f"Valid: {hmesh.valid(m)}")
    print(f"Closed: {hmesh.closed(m)}")

# Load mesh
m = hmesh.load(argv[1])
print_stats(m)

# Close holes and clean
hmesh.close_holes(m)
hmesh.triangulate(m)
hmesh.remove_caps(m)
hmesh.remove_needles(m)

# Print statistics after cleaning
print("After cleaning:")
print_stats(m)

# Smooth
hmesh.taubin_smooth(m,30)

# Save result
hmesh.save(argv[2], m)
```

### Mesh Simplification Pipeline

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl
from sys import argv

# Load high-resolution mesh
m = hmesh.load(argv[1])
original_faces = len(m.faces())

# Simplify to 5% of original faces
hmesh.quadric_simplify(m, keep_fraction=0.05)
simplified_faces = len(m.faces())

print(f"Reduced from {original_faces} to {simplified_faces} faces")

# Visualize result
viewer = gl.Viewer()
viewer.display(m, mode='w', bg_col=(1,1,1), smooth=False)
```

## Mesh Analysis

### Curvature Visualization

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl
from sys import argv

# Load mesh
m = hmesh.load(argv[1])

# Compute mean curvature at each vertex
curvatures = [m.mean_curvature(v) for v in m.vertices()]

# Display with color-coded curvature
viewer = gl.Viewer()
viewer.display(m, mode='s', data=curvatures)
```

### Topology Analysis

```python
import pygel3d.hmesh as hmesh
from sys import argv

def analyze_mesh(filename):
    m = hmesh.load(filename)
    
    # Topology
    components = hmesh.analyze_topology(m)
    for i,comp in enumerate(components):
        print(f"Component {i+1}:")
        print(f"- Vertices: {comp['V']}")
        print(f"- Edges: {comp['E']}")
        print(f"- Faces: {comp['F']}")
        print(f"- Number of boundary curves: {comp['b']}")
        print(f"- Genus: {comp['g']}")

    # Geometry
    bbox_min, bbox_max = hmesh.bbox(m)
    print(f"  Bounding box: {bbox_min} to {bbox_max}")
    
    sphere_center, sphere_radius = hmesh.bsphere(m)
    print(f"  Bounding sphere: center {sphere_center}, radius {sphere_radius}")
    
    print(f"  Total area: {hmesh.area(m):.4f}")
    if hmesh.closed(m):
        print(f"  Volume: {hmesh.volume(m):.4f}")
    
    # Edge statistics
    avg_len = hmesh.average_edge_length(m)
    med_len = hmesh.median_edge_length(m)
    print(f"  Average edge length: {avg_len:.4f}")
    print(f"  Median edge length: {med_len:.4f}")

# Analyze a mesh
analyze_mesh(argv[1])
```

## Graph Processing

### Skeleton Extraction and Visualization

```python
from sys import argv
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph
import pygel3d.gl_display as gl

# Load mesh
m = hmesh.load(argv[1])

# Extract graph
g = graph.from_mesh(m)
s = graph.MSLS_skeleton(g)


# Process skeleton
graph.prune(s)
l = s.average_edge_length()
graph.edge_contract(s, dist_thresh = l * 0.05)
graph.smooth(s, num_iter=10, alpha=0.1) 

# Convert to mesh for visualization
skeleton_mesh = hmesh.skeleton_to_feq(s, node_radii = 0.25 * l)

# Save
graph.save("skeleton.graph", s)
hmesh.save("skeleton.obj", skeleton_mesh)

# Visualize
viewer = gl.Viewer()
viewer.display(skeleton_mesh, s, mode='x')
```

## Spatial Queries

### Distance Field Computation

```python
from sys import argv
from pygel3d import hmesh
import numpy as np

# Load mesh
m = hmesh.load(argv[1])

# Create distance object
dist = hmesh.MeshDistance(m)

# Create sampling grid
grid_res = 50
bbox_min, bbox_max = hmesh.bbox(m)

# Expand bbox slightly
margin = 0.1
bbox_min = [x - margin for x in bbox_min]
bbox_max = [x + margin for x in bbox_max]

# Sample points
x = np.linspace(bbox_min[0], bbox_max[0], grid_res)
y = np.linspace(bbox_min[1], bbox_max[1], grid_res)
z = np.linspace(bbox_min[2], bbox_max[2], grid_res)

# Compute distance field
distances = np.zeros((grid_res, grid_res, grid_res))
for i, xi in enumerate(x):
    for j, yj in enumerate(y):
        for k, zk in enumerate(z):
            point = [xi, yj, zk]
            distances[i, j, k] = dist.signed_distance(point)

print(f"Distance field computed: {distances.shape}")
print(f"Min distance: {distances.min():.4f}")
print(f"Max distance: {distances.max():.4f}")

# Save distance field (e.g., as numpy array)
np.save("distance_field.npy", distances)
```

### Nearest Neighbor Search

```python
from pygel3d.spatial import I3DTree
import numpy as np
import time

# Generate random point cloud
n_points = 1000000
points = np.random.rand(n_points, 3) * 100

# Build kD-tree
print("Building tree...")
start = time.time()
tree = I3DTree()
for i, p in enumerate(points):
    tree.insert(p, i)
build_time = time.time() - start
print(f"Built tree with {n_points} points in {build_time:.2f}s")

# Query nearest neighbors
n_queries = 100000
query_points = np.random.rand(n_queries, 3) * 100

start = time.time()
for q in query_points:
    nearest_idx = tree.closest_point(q,1e20)
query_time = time.time() - start

print(f"Performed {n_queries} queries in {query_time:.2f}s")
print(f"Average query time: {query_time/n_queries*1000000:.2f}us")
```

## Advanced Workflows

### Mesh Reconstruction from Point Cloud

```python
import pygel3d.hmesh as hmesh
import numpy as np

# Load or generate point cloud
# points: Nx3 numpy array
# normals: Nx3 numpy array

# Example usage (with synthetic data)
theta = np.linspace(0, 2*np.pi, 100)[:-1]
phi = np.linspace(0.1, np.pi-0.1, 50)
theta, phi = np.meshgrid(theta, phi)

x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

points = np.stack([x.flatten(), y.flatten(), z.flatten()], axis=1)
normals = points / np.linalg.norm(points, axis=1, keepdims=True)

m = hmesh.rsr_recon(points, normals, max_normal_ang=10)
hmesh.close_holes(m, max_size=1000)
hmesh.triangulate(m)
for _ in range(5):
    hmesh.maximize_min_angle(m)
hmesh.save("reconstructed.obj", m)
```

## Jupyter Notebook Example

```python
# Cell 1: Setup
from pygel3d import hmesh, jupyter_display as jd

# Cell 2: Load and display original
m = hmesh.load("/path/to/model") # replace with your model path
print("Original mesh:")
jd.display(m)

# Cell 3: Smooth and display
hmesh.taubin_smooth(m)
print("After smoothing:")
jd.display(m)

# Cell 4: Simplify and display
hmesh.quadric_simplify(m, keep_fraction=0.5)
print("After simplification:")
jd.display(m)

# Cell 5: Save result
hmesh.save("result.obj", m)
print("Saved result.obj")
```

These examples demonstrate the most common PyGEL3D workflows. Modify them to suit your specific needs!
