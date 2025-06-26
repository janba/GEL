from pygel3d import hmesh, gl_display as gl

m = hmesh.load('../../../data/bunny.obj')

viewer = gl.Viewer()
viewer.display(m, smooth=False)

hmesh.close_holes(m, max_size=1000)
viewer.display(m, smooth=False)
hmesh.triangulate(m)
viewer.display(m, smooth=False)
hmesh.quadric_simplify(m, 0.1)
viewer.display(m, smooth=False)
hmesh.minimize_dihedral_angle(m)
viewer.display(m, smooth=False)
hmesh.maximize_min_angle(m, dihedral_thresh=0.98)
viewer.display(m, smooth=False)

