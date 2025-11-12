# Graph Module

::: pygel3d.graph

The `graph` module provides a 3D spatial graph data structure for representing curve skeletons and other graph-based geometric structures.

## Graph Class

The `Graph` class represents a spatial graph with 3D vertices connected by edges.

### Creating Graphs

```python
import pygel3d.graph as graph
import pygel3d.hmesh as hmesh

# Create empty graph
g = graph.Graph()

# Load from file
g = graph.load("skeleton.graph")

# Create from mesh
m = hmesh.load("model.obj")
g = graph.from_mesh(m)
```

## Graph I/O

### Loading and Saving
- `load(filename)` - Load graph from file
- `save(filename, graph)` - Save graph to file

## Graph Information

### Basic Queries
- `g.no_nodes()` - Number of nodes
- `g.no_edges()` - Number of edges
- `g.nodes()` - Get list of all node IDs
- `g.neighbors(node_id, mode)` - Get neighbors of a node

### Geometry
- `g.positions()` - Get all node positions as flat array
- `average_edge_length(g)` - Average edge length

## Graph Construction

### Adding Elements
- `g.add_node(position)` - Add a node at position [x, y, z]
- `g.connect_nodes(n0, n1)` - Connect two nodes with an edge

### Removing Elements
- `g.remove_node(node_id)` - Remove a node
- `g.disconnect_nodes(n0, n1)` - Remove edge between nodes

### Status
- `g.node_in_use(node_id)` - Check if node exists

## Graph Processing

### Cleaning
- `g.cleanup()` - Remove unused nodes
- `prune(g)` - Remove degree-1 nodes

### Smoothing
- `smooth(g, iterations, alpha)` - Smooth node positions

### Optimization
- `edge_contract(g, threshold)` - Contract short edges
- `saturate(g, hops, dist_fraction, radius)` - Add edges for connectivity

### Mesh Conversion
- `from_mesh(mesh)` - Extract graph from mesh
- `graph_to_mesh_cyl(g, mesh, fudge)` - Convert graph to cylindrical mesh

## Example Usage

### Graph Processing Pipeline

```python
import pygel3d.graph as graph
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Extract skeleton graph
g = graph.from_mesh(m)

# Process graph
graph.smooth(g, iter=10, alpha=0.5)
graph.prune(g)
graph.edge_contract(g, threshold=0.1)

# Save graph
graph.save("skeleton.graph", g)

# Convert to mesh for visualization
result = hmesh.Manifold()
graph.graph_to_mesh_cyl(g, result, fudge=0.5)
hmesh.save("skeleton_mesh.obj", result)
```

### Building a Custom Graph

```python
import pygel3d.graph as graph

# Create graph
g = graph.Graph()

# Add nodes
n0 = g.add_node([0, 0, 0])
n1 = g.add_node([1, 0, 0])
n2 = g.add_node([0.5, 1, 0])
n3 = g.add_node([0.5, 0.5, 1])

# Connect nodes
g.connect_nodes(n0, n1)
g.connect_nodes(n1, n2)
g.connect_nodes(n2, n3)
g.connect_nodes(n3, n0)

# Query
print(f"Nodes: {g.no_nodes()}")
print(f"Edges: {g.no_edges()}")
print(f"Neighbors of node 0: {g.neighbors(n0)}")
```
