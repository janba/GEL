# Usage: python mesh_simplification.py <input.obj>
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
