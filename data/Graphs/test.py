from pygel3d import hmesh, graph, gl_display as gl
s = graph.load("armadillo_symmetric.graph")
f = hmesh.skeleton_to_feq(s)
m = hmesh.load("../ReferenceMeshes/armadillo.obj")
hmesh.cc_split(f)
hmesh.cc_smooth(f)
hmesh.fit_mesh_to_ref(f, m)
v = gl.Viewer()
v.display(f)
