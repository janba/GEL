# Graph Processing Tutorial

This tutorial covers working with spatial graphs in PyGEL3D.

## Creating Graphs

### From Scratch

```python
import pygel3d.graph as graph

# Create empty graph
g = graph.Graph()

# Add nodes
n0 = g.add_node([0, 0, 0])
n1 = g.add_node([1, 0, 0])
n2 = g.add_node([0.5, 1, 0])

# Connect nodes
g.connect_nodes(n0, n1)
g.connect_nodes(n1, n2)
g.connect_nodes(n2, n0)

print(f"Created graph with {g.no_nodes()} nodes and {g.no_edges()} edges")
```

### From Mesh

```python
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph

# Load mesh
m = hmesh.load("model.obj")

# Extract skeleton graph
g = graph.from_mesh(m)
```

## Graph Queries

### Node Information

```python
# Get all node IDs
nodes = g.nodes()

# Get node positions
positions = g.positions()  # Flat array: [x0,y0,z0, x1,y1,z1, ...]

# Get neighbors
neighbors = g.neighbors(node_id, mode='n')

# Check if node exists
if g.node_in_use(node_id):
    print("Node exists")
```

### Graph Statistics

```python
print(f"Nodes: {g.no_nodes()}")
print(f"Edges: {g.no_edges()}")
print(f"Avg edge length: {graph.average_edge_length(g)}")
```

## Graph Processing

### Smoothing

```python
# Smooth graph positions
graph.smooth(g, iter=10, alpha=0.5)
```

### Pruning

```python
# Remove degree-1 nodes (leaf nodes)
graph.prune(g)
```

### Edge Contraction

```python
# Contract edges shorter than threshold
num_contracted = graph.edge_contract(g, threshold=0.1)
print(f"Contracted {num_contracted} edges")
```

## Graph to Mesh Conversion

### Cylindrical Mesh

```python
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph

# Load graph
g = graph.load("skeleton.graph")

# Convert to cylindrical mesh
m = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, m, fudge=0.5)

# Save result
hmesh.save("skeleton_mesh.obj", m)
```

## Complete Workflow

```python
import pygel3d.hmesh as hmesh
import pygel3d.graph as graph
import pygel3d.gl_display as gl

# 1. Load mesh
m = hmesh.load("model.obj")

# 2. Extract skeleton
g = graph.from_mesh(m)
print(f"Initial: {g.no_nodes()} nodes")

# 3. Process graph
graph.smooth(g, iter=10, alpha=0.5)
graph.prune(g)
contracted = graph.edge_contract(g, threshold=0.01)
print(f"After processing: {g.no_nodes()} nodes")

# 4. Save graph
graph.save("skeleton.graph", g)

# 5. Convert to mesh for visualization
skeleton_mesh = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, skeleton_mesh, fudge=0.5)

# 6. Visualize
viewer = gl.Viewer()
viewer.display(skeleton_mesh, mode='w')
```

See the [Graph API Reference](../api/graph.md) for complete documentation.
