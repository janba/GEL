from pygel3d import hmesh, graph, gl_display as gl
m = hmesh.load("as.obj")
s = graph.LS_skeleton(graph.from_mesh(m))
f = hmesh.skeleton_to_feq(s)
#hmesh.cc_split(f)
#hmesh.fit_mesh_to_ref(f, m)
v = gl.Viewer()
v.display(f)
