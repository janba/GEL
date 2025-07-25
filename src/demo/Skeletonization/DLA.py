from pygel3d import hmesh, gl_display as gl, graph
from numpy import where, array, linspace
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
rad_min = 0.1
rad_max = 1.0
g = graph.Graph()
g.add_node([.0,.0,.0])
for _ in range(10000):
    g.add_node(rad_min*random_point_in_unit_ball())
pos = g.positions()
fixed = array([True] + [False] * (len(pos) - 1))

# Diffusion Limited Aggregation (DLA) process:
# We will move the points randomly and connect them to the nearest fixed point.
# The process is repeated until all points are connected to a fixed point or we
# reach a maximum number of iterations.
for t in linspace(0, 1, 1500):
    pos = g.positions()
    rad = rad_min + (t**2) * (rad_max - rad_min)
    k = 0.01
    for i, x in enumerate(pos):
        if not fixed[i]:
            x += k * random_point_in_unit_ball()
            l = norm(x)
            if l > rad:
                x_norm = rad*x/l
                x -= 2*x_norm
    fixed_indices = where(fixed)[0]
    tree = cKDTree(pos[fixed])
    not_fixed_indices = where(~fixed)[0]
    nbors = tree.query_ball_point(pos[~fixed], 0.4*k, workers=-1)
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
graph.saturate(g, hops=3, rad=0.005)
# graph.smooth(g, iter=10, alpha=0.1)
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
graph.prune(s)
m = hmesh.graph_to_feq(s, node_radii=0.0002)
viewer2.display(m, smooth=False)
hmesh.cc_subdivide(m)
viewer2.display(m, mode='g')

# Save the mesh to an OBJ file.
hmesh.save("dla.obj", m)