# Expects model.obj in the current directory.
# Writes skeleton.graph and skeleton.obj.
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
