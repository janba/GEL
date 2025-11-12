# Visualization Tutorial

This tutorial covers visualization options in PyGEL3D.

## OpenGL Viewer

### Basic Display

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl

# Load mesh
m = hmesh.load("model.obj")

# Create viewer and display
viewer = gl.Viewer()
viewer.display(m)
```

### Rendering Modes

```python
# Wireframe
viewer.display(m, mode='w')

# Flat shading
viewer.display(m, mode='n')

# Smooth shading (glazed)
viewer.display(m, mode='g')

# Isophote lines
viewer.display(m, mode='i')

# Ghost (semi-transparent)
viewer.display(m, mode='x')
```

### Custom Background

```python
# White background
viewer.display(m, background=[1.0, 1.0, 1.0])

# Black background
viewer.display(m, background=[0.0, 0.0, 0.0])

# Gray background (default)
viewer.display(m, background=[0.3, 0.3, 0.3])
```

### Scalar Field Visualization

```python
import pygel3d.hmesh as hmesh
import pygel3d.gl_display as gl

# Load mesh
m = hmesh.load("model.obj")

# Compute scalar field (e.g., curvature)
curvatures = []
for v in m.vertices():
    curv = hmesh.mean_curvature(m, v)
    curvatures.append(abs(curv))

# Display with color-coding
viewer = gl.Viewer()
viewer.display(m, mode='s', data=curvatures)
```

### Vector Field Visualization

```python
# Compute vector field (e.g., normals)
normals = []
for v in m.vertices():
    normal = hmesh.vertex_normal(m, v)
    normals.extend(normal)  # Flatten

# Display as line field
viewer = gl.Viewer()
viewer.display(m, mode='l', data=normals)
```

## Jupyter Notebook Visualization

### Basic Display

```python
import pygel3d.hmesh as hmesh
import pygel3d.jupyter_display as jd

# Load and display
m = hmesh.load("model.obj")
jd.display(m)
```

### Custom Styling

```python
# With wireframe
jd.display(m, wireframe=True)

# Custom color
jd.display(m, color='coral')

# Custom size
jd.display(m, width=1000, height=800)

# Smooth or flat shading
jd.display(m, smooth=True)
```

### Side-by-Side Comparison

```python
from IPython.display import display, HTML

# Load two versions
m1 = hmesh.load("before.obj")
m2 = hmesh.load("after.obj")

# Display with labels
display(HTML("<h3>Before Processing</h3>"))
jd.display(m1, width=500, color='lightblue')

display(HTML("<h3>After Processing</h3>"))
jd.display(m2, width=500, color='lightgreen')
```

## Interactive Controls

### OpenGL Viewer Controls

- **Left Click + Drag**: Rotate camera
- **Right Click + Drag**: Zoom
- **Shift + Right Click + Drag**: Pan
- **ESC**: Exit viewer
- **Space**: Clear annotations
- **Ctrl + Click**: Add/remove annotation point

### Jupyter Widget Controls

- **Click + Drag**: Rotate
- **Scroll**: Zoom
- **Right Click + Drag**: Pan
- **Hover**: Show coordinates

## Multiple Views

### Multiple OpenGL Windows

```python
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

# Create multiple viewers
viewer1 = gl.Viewer()
viewer2 = gl.Viewer()

# Load meshes
m1 = hmesh.load("model1.obj")
m2 = hmesh.load("model2.obj")

# Display in different windows
viewer1.display(m1, mode='g')
viewer2.display(m2, mode='w')
```

## Tips and Best Practices

### For OpenGL Viewer

1. **Large Meshes**: Consider simplification for better performance
2. **Mode Selection**: Use 'w' (wireframe) to inspect topology
3. **Annotations**: Use Ctrl+Click to mark features
4. **Background**: Light backgrounds work well for screenshots

### For Jupyter Notebooks

1. **Export**: Notebooks export to HTML with interactive 3D
2. **Size**: Adjust width/height for better layout
3. **Colors**: Use contrasting colors for comparison
4. **Smooth Shading**: Usually looks better for organic shapes

### For Presentations

```python
# High-quality visualization setup
import pygel3d.gl_display as gl
import pygel3d.hmesh as hmesh

m = hmesh.load("model.obj")

# Clean mesh
hmesh.cc_smooth(m)

# Display with good settings
viewer = gl.Viewer()
viewer.display(
    m,
    mode='g',          # Smooth shading
    smooth=True,
    background=[1.0, 1.0, 1.0]  # White background
)
```

See the [GL Display API](../api/gl_display.md) and [Jupyter Display API](../api/jupyter_display.md) for complete documentation.
