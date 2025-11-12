# Examples

This page provides complete, runnable examples demonstrating common PyGEL3D workflows.

## Basic Mesh Processing

### Load, Process, and Save

```python
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("input.obj")

# Print statistics
print(f"Vertices: {m.no_vertices()}")
print(f"Faces: {m.no_faces()}")
print(f"Valid: {hmesh.valid(m)}")
print(f"Closed: {hmesh.closed(m)}")

# Clean and repair
hmesh.stitch_mesh(m, 1e-6)
hmesh.close_holes(m, 100)

# Smooth
hmesh.cc_smooth(m)

# Triangulate
hmesh.shortest_edge_triangulate(m)

# Save result
hmesh.save("output.obj", m)
```

### Mesh Simplification Pipeline

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl

# Load high-resolution mesh
m = hmesh.load("high_res_model.obj")
original_faces = m.no_faces()

# Simplify to 50% of original faces
hmesh.quadric_simplify(m, keep_fraction=0.5)
simplified_faces = m.no_faces()

print(f"Reduced from {original_faces} to {simplified_faces} faces")

# Visualize result
viewer = gl.Viewer()
viewer.display(m, mode='g')
```

## Mesh Analysis

### Curvature Visualization

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl

# Load mesh
m = hmesh.load("model.obj")

# Compute mean curvature at each vertex
curvatures = []
for v in m.vertices():
    curv = abs(hmesh.mean_curvature(m, v))
    curvatures.append(curv)

# Normalize for visualization
max_curv = max(curvatures)
if max_curv > 0:
    curvatures = [c/max_curv for c in curvatures]

# Display with color-coded curvature
viewer = gl.Viewer()
viewer.display(m, mode='s', data=curvatures)
```

### Topology Analysis

```python
import pygel3d.hmesh as hmesh

def analyze_mesh(filename):
    m = hmesh.load(filename)
    
    # Basic counts
    print(f"Mesh: {filename}")
    print(f"  Vertices: {m.no_vertices()}")
    print(f"  Faces: {m.no_faces()}")
    print(f"  Halfedges: {m.no_halfedges()}")
    
    # Topology
    print(f"  Valid: {hmesh.valid(m)}")
    print(f"  Closed: {hmesh.closed(m)}")
    print(f"  Boundaries: {hmesh.count_boundary_curves(m)}")
    
    # Geometry
    bbox_min, bbox_max = hmesh.bbox(m)
    print(f"  Bounding box: {bbox_min} to {bbox_max}")
    
    sphere_center, sphere_radius = hmesh.bsphere(m)
    print(f"  Bounding sphere: center {sphere_center}, radius {sphere_radius}")
    
    print(f"  Total area: {hmesh.total_area(m):.4f}")
    if hmesh.closed(m):
        print(f"  Volume: {hmesh.volume(m):.4f}")
    
    # Edge statistics
    avg_len = hmesh.average_edge_length(m)
    med_len = hmesh.median_edge_length(m)
    print(f"  Average edge length: {avg_len:.4f}")
    print(f"  Median edge length: {med_len:.4f}")

# Analyze a mesh
analyze_mesh("bunny.obj")
```

## Graph Processing

### Skeleton Extraction and Visualization

```python
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph
import pygel3d.gl_display as gl

# Load mesh
m = hmesh.load("model.obj")

# Extract skeleton
g = graph.from_mesh(m)

print(f"Skeleton: {g.no_nodes()} nodes, {g.no_edges()} edges")

# Process skeleton
graph.smooth(g, iter=10, alpha=0.5)
graph.prune(g)
graph.edge_contract(g, threshold=0.01)

print(f"After processing: {g.no_nodes()} nodes, {g.no_edges()} edges")

# Convert to mesh for visualization
skeleton_mesh = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, skeleton_mesh, fudge=0.5)

# Save
graph.save("skeleton.graph", g)
hmesh.save("skeleton.obj", skeleton_mesh)

# Visualize
viewer = gl.Viewer()
viewer.display(skeleton_mesh, mode='w')
```

## Spatial Queries

### Distance Field Computation

```python
import pygel3d.hmesh as hmesh
from pygel3d import MeshDistance
import numpy as np

# Load mesh
m = hmesh.load("bunny.obj")

# Create distance object
dist = MeshDistance(m)

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
from pygel3d import I3DTree
import numpy as np
import time

# Generate random point cloud
n_points = 100000
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
n_queries = 1000
query_points = np.random.rand(n_queries, 3) * 100

start = time.time()
for q in query_points:
    nearest_idx = tree.closest_point(q)
query_time = time.time() - start

print(f"Performed {n_queries} queries in {query_time:.2f}s")
print(f"Average query time: {query_time/n_queries*1000:.2f}ms")
```

## Advanced Workflows

### Mesh Reconstruction from Point Cloud

```python
import pygel3d.hmesh as hmesh
import numpy as np

# Load or generate point cloud
# points: Nx3 numpy array
# normals: Nx3 numpy array

def reconstruct_mesh(points, normals):
    # Flatten arrays
    pts_flat = points.flatten().tolist()
    nrm_flat = normals.flatten().tolist()
    
    # Create mesh
    m = hmesh.Manifold()
    
    # Reconstruct using RSR
    # Parameters need to be tuned for specific data
    hmesh.rsr_recon(
        m,
        pts_flat,
        nrm_flat,
        len(points),
        len(normals),
        isEuclidean=True,
        genus=0,
        k=4,
        r=2,
        theta=30,
        n=10
    )
    
    return m

# Example usage (with synthetic data)
theta = np.linspace(0, 2*np.pi, 100)
phi = np.linspace(0, np.pi, 50)
theta, phi = np.meshgrid(theta, phi)

x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

points = np.stack([x.flatten(), y.flatten(), z.flatten()], axis=1)
normals = points / np.linalg.norm(points, axis=1, keepdims=True)

m = reconstruct_mesh(points, normals)
hmesh.save("reconstructed.obj", m)
```

### Batch Processing

```python
import pygel3d.hmesh as hmesh
import os
from pathlib import Path

def process_mesh_directory(input_dir, output_dir, keep_fraction=0.5):
    """Process all meshes in a directory."""
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Find all mesh files
    mesh_files = list(input_path.glob("*.obj"))
    
    print(f"Found {len(mesh_files)} meshes to process")
    
    for i, mesh_file in enumerate(mesh_files):
        print(f"Processing {i+1}/{len(mesh_files)}: {mesh_file.name}")
        
        try:
            # Load
            m = hmesh.load(str(mesh_file))
            
            # Process
            hmesh.stitch_mesh(m, 1e-6)
            hmesh.close_holes(m, 100)
            hmesh.cc_smooth(m)
            hmesh.shortest_edge_triangulate(m)
            hmesh.quadric_simplify(m, keep_fraction)
            
            # Save
            output_file = output_path / mesh_file.name
            hmesh.save(str(output_file), m)
            
            print(f"  Saved to {output_file}")
            
        except Exception as e:
            print(f"  Error processing {mesh_file.name}: {e}")

# Process directory
process_mesh_directory("input_meshes", "output_meshes", keep_fraction=0.5)
```

## Jupyter Notebook Example

```python
# Cell 1: Setup
import pygel3d.hmesh as hmesh
import pygel3d.jupyter_display as jd

# Cell 2: Load and display original
m = hmesh.load("model.obj")
print("Original mesh:")
jd.display(m, color='lightblue')

# Cell 3: Smooth and display
hmesh.cc_smooth(m)
print("After smoothing:")
jd.display(m, color='lightgreen')

# Cell 4: Simplify and display
hmesh.quadric_simplify(m, keep_fraction=0.5)
print("After simplification:")
jd.display(m, color='coral')

# Cell 5: Save result
hmesh.save("result.obj", m)
print("Saved result.obj")
```

These examples demonstrate the most common PyGEL3D workflows. Modify them to suit your specific needs!
