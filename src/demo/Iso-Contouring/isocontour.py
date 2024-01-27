from pygel3d import hmesh, graph, gl_display as gl

import numpy as np
from math import cos, sin
array = np.array

v = gl.Viewer()

dim = (40,40,40)
vol = np.zeros(dim)
for idx in np.ndindex(dim):
    x,y,z = array(idx) * 0.25
    vol[idx] = sin(x)*cos(y)+sin(y)*cos(z)+sin(z)*cos(x)

m = hmesh.volumetric_isocontour(vol, high_is_inside=True)
v.display(m, mode='w', smooth=True)

m = hmesh.volumetric_isocontour(vol, high_is_inside=True, dual_connectivity=True)
v.display(m, mode='w', smooth=True)


dim = (16,16,16)
vol = np.zeros(dim)
n = array((.5,.3,.1))
p0 = array((6,6,6))
for idx in np.ndindex(dim):
    p = array(idx)
    vol[idx] = n @ (p - p0)

m = hmesh.volumetric_isocontour(vol, make_triangles=False, high_is_inside=False)
v.display(m,smooth=False)

m_in = hmesh.load("../../../data/bunny.obj")
v.display(m_in, reset_view=True)
D = hmesh.MeshDistance(m_in)
plo, phi = hmesh.bbox(m_in)

dim = tuple(int(x) for x in np.ceil((phi-plo) / 0.0025))
spacing = (phi-plo)/(array(dim) - 1)
vol = np.array([D.signed_distance(array(idx)*spacing + plo) for idx in np.ndindex(dim)])
vol = vol.reshape(dim)

m = hmesh.volumetric_isocontour(vol, plo, phi, high_is_inside=False)
d = np.clip(D.signed_distance(m.positions()), -0.0005, 0.0005)
v.display(m, mode='w')
v.display(m, mode='s', data=d)
