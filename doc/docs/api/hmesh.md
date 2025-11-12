# HMesh Module

::: pygel3d.hmesh

The `hmesh` module provides the core halfedge mesh data structure and associated operations for polygonal mesh processing.

## Manifold Class

The `Manifold` class represents a polygonal mesh using the halfedge data structure, which enables efficient traversal and manipulation of mesh topology.

### Creating Meshes

You can create meshes in several ways:

```python
import pygel3d.hmesh as hmesh

# Create an empty mesh
m = hmesh.Manifold()

# Load from file
m = hmesh.load("model.obj")
m = hmesh.obj_load("model.obj")
m = hmesh.off_load("model.off")
m = hmesh.ply_load("model.ply")
m = hmesh.x3d_load("model.x3d")
```

## Mesh I/O Functions

### Loading Meshes
- `load(filename, mesh)` - Load mesh from file (auto-detect format)
- `obj_load(filename, mesh)` - Load Wavefront OBJ file
- `off_load(filename, mesh)` - Load Object File Format
- `ply_load(filename, mesh)` - Load PLY file
- `x3d_load(filename, mesh)` - Load X3D file

### Saving Meshes
- `obj_save(filename, mesh)` - Save as Wavefront OBJ
- `off_save(filename, mesh)` - Save as Object File Format
- `x3d_save(filename, mesh)` - Save as X3D

## Mesh Information

### Basic Queries
- `valid(mesh)` - Check if mesh is valid
- `closed(mesh)` - Check if mesh is closed (no boundary)
- `bbox(mesh)` - Get bounding box (returns min, max)
- `bsphere(mesh)` - Get bounding sphere (returns center, radius)

### Counts
- `mesh.no_vertices()` - Number of vertices
- `mesh.no_faces()` - Number of faces
- `mesh.no_halfedges()` - Number of halfedges
- `count_boundary_curves(mesh)` - Count boundary loops

## Mesh Processing

### Cleaning and Repair
- `stitch_mesh(mesh, threshold)` - Merge nearby vertices
- `close_holes(mesh, max_size)` - Fill holes up to max_size edges
- `remove_caps(mesh, threshold)` - Remove cap-like features
- `remove_needles(mesh, threshold, avg_pos)` - Remove needle-like features
- `merge_coincident_boundary_vertices(mesh, threshold)` - Merge boundary vertices

### Smoothing
- `cc_smooth(mesh)` - Catmull-Clark smoothing
- `loop_smooth(mesh)` - Loop smoothing
- `taubin_smooth(mesh, iterations)` - Taubin smoothing
- `laplacian_smooth(mesh, weight, iterations)` - Laplacian smoothing
- `anisotropic_smooth(mesh, sharpness, iterations)` - Anisotropic smoothing

### Subdivision
- `cc_split(mesh)` - Catmull-Clark subdivision
- `loop_split(mesh)` - Loop subdivision (for triangles)
- `root3_subdivide(mesh)` - Root-3 subdivision
- `butterfly_subdivide(mesh)` - Butterfly subdivision

### Simplification
- `quadric_simplify(mesh, keep_fraction, singular_threshold, error_threshold)` - Quadric error metric simplification

### Refinement
- `refine_edges(mesh, threshold)` - Refine long edges

### Triangulation
- `shortest_edge_triangulate(mesh)` - Triangulate using shortest diagonal
- `ear_clip_triangulate(mesh)` - Triangulate using ear clipping

### Optimization
- `minimize_curvature(mesh, anneal)` - Minimize mesh curvature
- `minimize_dihedral_angle(mesh, max_iter, anneal, alpha, gamma)` - Minimize dihedral angles
- `maximize_min_angle(mesh, threshold, anneal)` - Maximize minimum angle
- `optimize_valency(mesh, anneal)` - Optimize vertex valency
- `randomize_mesh(mesh, max_iter)` - Random edge flips

### Topology Operations
- `flip_orientation(mesh)` - Reverse face orientation
- `cc_split(mesh)` - Catmull-Clark split

## Mesh Measurements

### Geometric Measurements
- `area(mesh, face_id)` - Area of a face
- `perimeter(mesh, face_id)` - Perimeter of a face
- `length(mesh, halfedge_id)` - Length of an edge
- `total_area(mesh)` - Total surface area
- `volume(mesh)` - Mesh volume

### Vertex Measurements
- `valency(mesh, vertex_id)` - Vertex valency (degree)
- `one_ring_area(mesh, vertex_id)` - Area of one-ring neighborhood
- `mixed_area(mesh, vertex_id)` - Mixed Voronoi/barycentric area

### Curvature
- `gaussian_curvature(mesh, vertex_id)` - Gaussian curvature at vertex
- `mean_curvature(mesh, vertex_id)` - Mean curvature at vertex
- `principal_curvatures(mesh, vertex_id)` - Principal curvatures

### Normals
- `vertex_normal(mesh, vertex_id)` - Vertex normal
- `face_normal(mesh, face_id)` - Face normal

## Mesh Traversal

### Vertex Operations
- `mesh.vertices()` - Iterator over all vertices
- `mesh.circulate_vertex(vertex_id, mode)` - Circulate around vertex

### Face Operations
- `mesh.faces()` - Iterator over all faces
- `mesh.circulate_face(face_id, mode)` - Circulate around face
- `no_edges(mesh, face_id)` - Number of edges in face
- `centre(mesh, face_id)` - Face center

### Halfedge Operations
- `mesh.halfedges()` - Iterator over all halfedges

## Mesh Editing

### Adding Elements
- `mesh.add_vertex(position)` - Add a vertex
- `mesh.add_face(positions)` - Add a face

### Removing Elements
- `mesh.remove_vertex(vertex_id)` - Remove a vertex
- `mesh.remove_face(face_id)` - Remove a face
- `mesh.remove_edge(halfedge_id)` - Remove an edge

### Modifying Topology
- `mesh.flip_edge(halfedge_id)` - Flip an edge
- `mesh.collapse_edge(halfedge_id, avg_vertices)` - Collapse an edge
- `mesh.split_edge(halfedge_id)` - Split an edge
- `mesh.split_face_by_edge(face_id, v0, v1)` - Split face with new edge
- `mesh.split_face_by_vertex(face_id)` - Split face from center
- `mesh.merge_faces(face_id, halfedge_id)` - Merge two faces
- `mesh.stitch_boundary_edges(h0, h1)` - Stitch two boundary edges
- `mesh.close_hole(halfedge_id)` - Close a boundary loop

### Status Queries
- `mesh.vertex_in_use(vertex_id)` - Check if vertex exists
- `mesh.face_in_use(face_id)` - Check if face exists
- `mesh.halfedge_in_use(halfedge_id)` - Check if halfedge exists

### Boundary Queries
- `is_vertex_at_boundary(mesh, vertex_id)` - Check if vertex is on boundary
- `is_halfedge_at_boundary(mesh, halfedge_id)` - Check if halfedge is on boundary
- `boundary_edge(mesh, vertex_id, halfedge_id)` - Check if edge is on boundary

### Connectivity
- `connected(mesh, v0, v1)` - Check if two vertices are connected

## Walker Functions

Walker functions provide low-level halfedge traversal:

- `mesh.walker.next_halfedge(h)` - Get next halfedge in face
- `mesh.walker.prev_halfedge(h)` - Get previous halfedge in face
- `mesh.walker.opposite_halfedge(h)` - Get opposite halfedge
- `mesh.walker.incident_face(h)` - Get incident face
- `mesh.walker.incident_vertex(h)` - Get incident vertex

## Advanced Operations

### Volumetric Operations
- `volumetric_isocontour(mesh, x_dim, y_dim, z_dim, data, pmin, pmax, tau, make_triangles, high_is_inside, dual_connectivity)` - Extract isosurface from volume

### Registration
- `non_rigid_registration(mesh, reference_mesh)` - Non-rigid registration
- `stable_marriage_registration(mesh, reference_mesh)` - Stable marriage registration

### Reconstruction
- `rsr_recon(mesh, vertices, normals, v_num, n_num, isEuclidean, genus, k, r, theta, n)` - Rotation system reconstruction

### Face Operations
- `extrude_faces(mesh, face_list, output_face_list)` - Extrude faces
- `kill_face_loop(mesh)` - Remove face loop
- `kill_degenerate_face_loops(mesh, threshold)` - Remove degenerate loops

### Connected Components
- `connected_components(mesh)` - Split into connected components

## MeshDistance Class

The `MeshDistance` class provides efficient distance queries to triangle meshes.

```python
from pygel3d import MeshDistance
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Create distance object
dist = MeshDistance(m)

# Query signed distance
point = [0, 0, 0]
distance = dist.signed_distance(point)

# Query unsigned distance
distance = dist.distance(point)
```

## Example Usage

### Complete Mesh Processing Pipeline

```python
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("input.obj")

# Clean mesh
hmesh.stitch_mesh(m, 1e-6)
hmesh.close_holes(m, 100)
hmesh.remove_caps(m, 0.1)

# Smooth
hmesh.cc_smooth(m)

# Triangulate
hmesh.shortest_edge_triangulate(m)

# Simplify
hmesh.quadric_simplify(m, keep_fraction=0.5)

# Optimize
hmesh.minimize_curvature(m, anneal=True)

# Save
hmesh.obj_save("output.obj", m)
```
