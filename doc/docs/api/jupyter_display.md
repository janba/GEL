# Jupyter Display Module

::: pygel3d.jupyter_display

The `jupyter_display` module provides visualization tools for Jupyter notebooks using Plotly, enabling interactive 3D graphics that can be exported to HTML.

## Display Function

The main function for displaying geometry in Jupyter notebooks:

```python
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Display in notebook
jd.display(m, smooth=True, wireframe=False)
```

## Function Parameters

### display(obj, **kwargs)

- `obj`: Manifold or Graph object to display
- `smooth`: Use smooth shading (default: True)
- `wireframe`: Show wireframe (default: False)
- `color`: Mesh color as string or RGB list (default: 'lightblue')
- `width`: Figure width in pixels (default: 800)
- `height`: Figure height in pixels (default: 600)

## Features

### Interactive 3D Widgets
- Rotate, pan, and zoom with mouse
- Hover to see coordinates
- Camera controls in toolbar
- Full Plotly interactivity

### HTML Export
- Notebooks can be exported to HTML
- 3D visualizations remain interactive
- Perfect for assignments and presentations
- Works in nbviewer and GitHub

### Google Colab Support
- Fully compatible with Google Colab
- No additional setup required
- Same functionality as local Jupyter

## Example Usage

### Basic Display

```python
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("bunny.obj")

# Display
jd.display(m)
```

### Custom Styling

```python
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh

# Load mesh
m = hmesh.load("model.obj")

# Display with custom appearance
jd.display(m, 
           smooth=True,
           wireframe=True,
           color='coral',
           width=1000,
           height=800)
```

### Multiple Objects

```python
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh

# Load multiple meshes
m1 = hmesh.load("model1.obj")
m2 = hmesh.load("model2.obj")

# Display separately
print("Model 1:")
jd.display(m1, color='red')

print("Model 2:")
jd.display(m2, color='blue')
```

### Color by Attribute

```python
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh
import numpy as np

# Load mesh
m = hmesh.load("model.obj")

# Compute per-vertex attribute (e.g., height)
positions = m.positions()
z_coords = [positions[3*i+2] for i in range(m.no_vertices())]

# Display with color mapping
# Note: For custom attribute coloring, you may need to create
# a custom Plotly figure or use the hmesh scalar field functions
jd.display(m)
```

### Displaying Graphs

```python
import pygel3d.jupyter_display as jd
import pygel3d.graph as graph

# Load graph
g = graph.load("skeleton.graph")

# Display
jd.display(g, color='green')
```

## Jupyter Notebook Workflow

### Setup Cell

```python
# Install if needed
!pip install PyGEL3D plotly

# Import modules
import pygel3d.hmesh as hmesh
import pygel3d.jupyter_display as jd
```

### Processing and Visualization

```python
# Load mesh
m = hmesh.load("input.obj")

# Show original
print("Original mesh:")
jd.display(m)

# Process
hmesh.cc_smooth(m)
hmesh.triangulate(m)

# Show result
print("Processed mesh:")
jd.display(m, color='lightgreen')
```

### Comparison View

```python
# Load two meshes
m1 = hmesh.load("before.obj")
m2 = hmesh.load("after.obj")

# Display side by side
from IPython.display import display, HTML

display(HTML("<h3>Before</h3>"))
jd.display(m1, width=400)

display(HTML("<h3>After</h3>"))
jd.display(m2, width=400)
```

## Google Colab Setup

```python
# First cell: Install dependencies
!apt-get install libglu1 libgl1
!pip install PyGEL3D plotly

# Import modules
import pygel3d.hmesh as hmesh
import pygel3d.jupyter_display as jd

# Upload file (if needed)
from google.colab import files
uploaded = files.upload()

# Load and display
m = hmesh.load(list(uploaded.keys())[0])
jd.display(m)
```

## Tips and Best Practices

### Performance
- Large meshes may be slow in the browser
- Consider simplification for complex models
- Use decimation for very large meshes

### Visualization
- Smooth shading is better for organic models
- Wireframe helps see topology
- Light colors work well on default backgrounds

### Notebooks
- Add markdown cells to explain each step
- Use clear section headings
- Include parameter descriptions
- Export to HTML for sharing

### Export Quality
- Ensure cells are executed before export
- Test exported HTML in browser
- Check that 3D widgets are interactive
- Use "Trust Notebook" if needed

## Plotly Integration

The module uses Plotly, so you can access advanced Plotly features:

```python
import plotly.graph_objects as go
import pygel3d.jupyter_display as jd
import pygel3d.hmesh as hmesh

# For advanced customization, you may need to work
# directly with Plotly's mesh3d objects
# See Plotly documentation for details
```

## Troubleshooting

### Widget Not Displaying
- Ensure Plotly is installed: `pip install plotly`
- Restart kernel and re-run cells
- Check browser console for errors

### Export Issues
- Use "File > Download as > HTML" in Jupyter
- Ensure all cells are executed
- Check that notebook is trusted

### Performance Issues
- Simplify large meshes before display
- Reduce figure size (width/height)
- Close unused notebooks
