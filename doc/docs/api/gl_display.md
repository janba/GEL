# GL Display Module

::: pygel3d.gl_display

The `gl_display` module provides OpenGL-based interactive 3D visualization for meshes and graphs.

## Viewer Class

The `Viewer` class creates an OpenGL window for displaying and interacting with 3D geometry.

### Creating a Viewer

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Create viewer
viewer = gl.Viewer()

# Load and display mesh
m = hmesh.load("model.obj")
viewer.display(m)
```

## Display Methods

### Main Display
- `viewer.display(mesh, mode, smooth, background, data)` - Display a mesh or graph

### Parameters
- `mesh`: Manifold or Graph object to display
- `mode`: Rendering mode (see below)
- `smooth`: Enable smooth shading (default: True)
- `background`: Background color as [r, g, b] (default: [0.3, 0.3, 0.3])
- `data`: Optional attribute data for scalar/vector field visualization

## Rendering Modes

The `mode` parameter controls how the geometry is rendered:

- `'w'`: **Wireframe** - Show edges only
- `'n'`: **Normal** - Flat shading with face normals
- `'g'`: **Glazed** - Smooth shading (default)
- `'i'`: **Isophote** - Isophote lines for curvature analysis
- `'l'`: **Line Field** - Display vector field as lines (requires data)
- `'s'`: **Scalar Field** - Color-coded scalar values (requires data)
- `'x'`: **Ghost** - Semi-transparent rendering

## Interactive Controls

### Mouse Controls
- **Left Mouse**: Rotate camera
- **Right Mouse**: Zoom in/out
- **Shift + Right Mouse**: Pan camera

### Keyboard Controls
- **ESC**: Exit viewer
- **Space**: Clear annotations

## Annotation

The viewer supports interactive point annotation:

- **Ctrl + Click**: Add/remove annotation point
- Annotation points are displayed as colored spheres
- Useful for marking features or measurements

## Multiple Viewers

You can create multiple viewers to display different objects:

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Create two viewers
viewer1 = gl.Viewer()
viewer2 = gl.Viewer()

# Display different meshes
m1 = hmesh.load("model1.obj")
m2 = hmesh.load("model2.obj")

viewer1.display(m1, mode='g')
viewer2.display(m2, mode='w')
```

## Example Usage

### Basic Visualization

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("bunny.obj")

# Create viewer and display
viewer = gl.Viewer()
viewer.display(m, mode='g', smooth=True)
```

### Custom Background

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Display with white background
viewer = gl.Viewer()
viewer.display(m, 
               mode='g', 
               smooth=True, 
               background=[1.0, 1.0, 1.0])
```

### Scalar Field Visualization

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh
import numpy as np

# Load mesh
m = hmesh.load("model.obj")

# Compute scalar field (e.g., mean curvature)
n_verts = m.no_vertices()
curvatures = []
for v in m.vertices():
    curv = hmesh.mean_curvature(m, v)
    curvatures.append(curv)

# Display with color-coded curvature
viewer = gl.Viewer()
viewer.display(m, mode='s', data=curvatures)
```

### Vector Field Visualization

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Compute vector field (e.g., normals)
n_verts = m.no_vertices()
normals = []
for v in m.vertices():
    normal = hmesh.vertex_normal(m, v)
    normals.extend(normal)  # Flatten to [x0,y0,z0,x1,y1,z1,...]

# Display with line field
viewer = gl.Viewer()
viewer.display(m, mode='l', data=normals)
```

### Wireframe Overlay

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Display as wireframe for topology inspection
viewer = gl.Viewer()
viewer.display(m, mode='w', background=[1.0, 1.0, 1.0])
```

### Displaying Graphs

```python
import pygel3d.gl_display as gl
import pygel3d.graph as graph

# Load or create graph
g = graph.load("skeleton.graph")

# Display graph
viewer = gl.Viewer()
viewer.display(g)
```

## Tips

- **Performance**: Large meshes may be slow to render; consider simplification
- **Smooth Shading**: Better for organic shapes; use flat for architectural models
- **Background**: Light backgrounds work well for presentations
- **Mode Selection**: Try different modes to highlight different features
- **Annotations**: Use for measurements or marking regions of interest
