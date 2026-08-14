# Usage: python curvature_visualization.py <input.obj>
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
