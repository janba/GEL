"""This script performs very basic test of the PyGEL API."""
from pygel3d import gl_display as gl, hmesh
from random import uniform

m = hmesh.load("../../../data/tetra.obj")

v = gl.Viewer()
v.display(m, mode='w')
for i in m.vertices():
    m.positions()[i] += 0.5 * m.vertex_normal(i)
v.display(m, mode='w')
hmesh.bsphere(m)
hmesh.bbox(m)


pts = [[-1.0, -1.0, -1.0], 
       [1.0, -1.0, -1.0], 
       [1.0, 1.0, -1.0], 
       [-1.0, 1.0, -1.0]]

tri = [[0, 1, 2], 
       [0, 2, 3]]

m = hmesh.Manifold.from_triangles(pts, tri)
v.display(m, mode='w')

pts = [ [uniform(0,1), uniform(0,1), 0.0] for i in range(100)]

m = hmesh.Manifold.from_points(pts)
v.display(m, mode='w')
