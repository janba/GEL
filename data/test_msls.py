from pygel3d import hmesh, gl_display as gl, graph
m = hmesh.load("Armadillo.ply")
g = graph.from_mesh(m)
s = graph.MSLS_skeleton(g, grow_thresh=512)
v = gl.Viewer()
v.display(m, s, mode="x")
