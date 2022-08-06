from pygel3d import hmesh, graph, gl_display as gl

s = graph.load('../../../data/Graphs/bunny.graph')

m_skel = hmesh.skeleton_to_feq(s)

print('Inverse Skeletonized')

ref_mesh = hmesh.load('../../../data/ReferenceMeshes/bunny.obj')

hmesh.cc_split(m_skel)
hmesh.cc_smooth(m_skel)

print('Subdivided, now fitting.')

fit_mesh = m_skel

fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, 50)
 
viewer = gl.Viewer()

viewer.display(fit_mesh)


