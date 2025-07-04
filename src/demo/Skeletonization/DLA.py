from pygel3d import hmesh, gl_display as gl, graph
from numpy import where, array
from numpy.linalg import norm
from scipy.spatial import cKDTree
from numpy.random import uniform

def random_point_in_unit_ball():
    """Generate a random point in the unit ball."""
    while True:
        x = uniform(-1.0, 1.0, size=3)
        if norm(x) < 1.0:
            return x

# Create a viewer.
viewer = gl.Viewer()

# Create a graph and add a central node at the origin.
# Then add 5000 random points in the unit ball around the origin.
# The points are added as nodes to the graph.
g = graph.Graph()
g.add_node([.0,.0,.0])
for _ in range(10000):
    g.add_node(random_point_in_unit_ball())
pos = g.positions()
fixed = array([True] + [False] * (len(pos) - 1))

# Diffusion Limited Aggregation (DLA) process:
# We will move the points randomly and connect them to the nearest fixed point.
# The process is repeated until all points are connected to a fixed point or we
# reach a maximum number of iterations.
for _ in range(1000):
    for i, x in enumerate(pos):
        if not fixed[i]:
            x += 0.025 * random_point_in_unit_ball()
            l = norm(x)
            if l > 1.0:
                x_norm = x/l
                x -= 2*x_norm
    fixed_indices = where(fixed)[0]
    tree = cKDTree(pos[fixed])
    not_fixed_indices = where(~fixed)[0]
    nbors = tree.query_ball_point(pos[~fixed], 0.035, workers=-1)
    for i, nbor_list in zip(not_fixed_indices, nbors):
        if len(nbor_list) > 0:
            fixed[i] = True
            for j in nbor_list:
                g.connect_nodes(i,fixed_indices[j])
    viewer.display(g, once=True)

# Cleanup: remove nodes that are not fixed and add some edges to the graph (saturate).
# We show the resulting graph in the viewer and compute the MSLS skeleton.
for i in g.nodes():
    if not fixed[i]:
        g.remove_node(i)
graph.saturate(g, hops=3, rad=0.05)
viewer.display(g)
s = graph.MSLS_skeleton(g)

# Display the skeleton in a new viewer. We clone the controller from the first viewer
# so that the view is the same.
viewer2 = gl.Viewer()
viewer2.clone_controller(viewer)
viewer2.display(s)

# Convert the skeleton to a manifold mesh and display it. node_radii is the thickness
# of the skeleton edges when converted to a mesh.
# The mesh is displayed first in an opaque fashion and the a second time in x-ray
# mode to show the skeleton edges inside the mesh.
m = hmesh.graph_to_feq(s, node_radii=0.005)
viewer2.display(m, smooth=False)
viewer2.display(m, s, mode='x', reset_view=True)

# Save the mesh to an OBJ file.
hmesh.save("dla.obj", m)