#!/opt/local/bin/python
from pygel3d import hmesh, gl_display as gl
from numpy.linalg import norm
from numpy import array, float64

def extrude(m, faces):
    n = array([0.0,0.0,0.0], dtype=float64)
    a = 0
    for f in faces:
        n += m.face_normal(f)
        a += m.area(f)
    a /= 4
    n /= norm(n)
    fset_out = hmesh.extrude_faces(m,faces)
    vset = set()
    pos = m.positions()
    for f in faces:
        verts =  m.circulate_face(f)
        for v in verts:
            vset.add(v)
    avg_pos = array([0.0,0,0])
    for v in vset:
        pos[v] += (a+0.15)*n
        avg_pos += pos[v]
    avg_pos /= len(vset)
    for v in vset:
        pos[v] = 0.25 * pos[v] + 0.75 * avg_pos
    return fset_out

viewer = gl.Viewer()
m = hmesh.load("../../../data/cube.obj")

fset_in = [0,3]

for _ in range(3):
    fset_out = set()
    for f in fset_in:
        fset_out |= extrude(m, [f])
    fset_in = list(fset_out)

viewer.display(m)



# m = hmesh.load("../../../data/arma_quads.obj")
# viewer.display(m)
# while True:
#     hmesh.kill_face_loop(m)
#     viewer.display(m, once=True)

