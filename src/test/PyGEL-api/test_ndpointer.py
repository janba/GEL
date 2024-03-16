from pygel3d import gl_display as gl, hmesh
m = hmesh.load("tetra.obj")
v = gl.Viewer()
v.display(m, mode='g')
for i in m.vertices():
    m.positions()[i] += 0.5 * m.vertex_normal(i)
v.display(m, mode='g')
hmesh.bsphere(m)
hmesh.bbox(m)
print(m.centre(0))
