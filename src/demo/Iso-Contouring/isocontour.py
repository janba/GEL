from math import cos, sin
from numpy import array, zeros, ndindex, ceil, clip
import numpy as np
from pygel3d import hmesh, graph, gl_display as gl

v = gl.Viewer()

# First test, create a volume containing a gyroid and display its isocontour
dim = (40,50,60)
vol = np.zeros(dim)
for idx in np.ndindex(dim):
    x,y,z = [ i * 0.25 for i in idx ]
    vol[idx] = sin(x)*cos(y)+sin(y)*cos(z)+sin(z)*cos(x)
m = hmesh.volumetric_isocontour(vol, high_is_inside=True)
v.display(m, mode='g', smooth=False)

# Second test, create a volume containing a gyroid but cut out some slices by setting
# the voxels to NaN and display the isocontour. Surfaces are now open.
vol[:,20:30,:] = float('nan')
vol[15:25,:,:] = float('nan')
vol[:,:,25:35] = float('nan')
m = hmesh.volumetric_isocontour(vol, high_is_inside=True)
v.display(m, mode='g', smooth=False)

# Third test, create a volume containing a half space delimited by the volume and 
# display the isocontour. This time the surface is closed.
dim = (16,16,16)
vol = np.zeros(dim)
n = array((.5,.3,.1))
p0 = array((6,6,6))
for idx in np.ndindex(dim):
    p = array(idx)
    vol[idx] = n @ (p - p0)
m = hmesh.volumetric_isocontour(vol, make_triangles=False, high_is_inside=False)
v.display(m,smooth=False,reset_view=True)

# Fourth test, create a volume from a signed distance field of the bunny and 
# display the mesh contoured using dual contouring. Show also the distance field
# as a scalar field on the mesh.
m_in = hmesh.load("../../../data/bunny.obj")
D = hmesh.MeshDistance(m_in)
plo, phi = hmesh.bbox(m_in)
dim = tuple(int(x) for x in np.ceil((phi-plo) / 0.0025))
spacing = (phi-plo)/(array(dim) - 1)
vol = D.signed_distance( [ array(idx)*spacing + plo for idx in np.ndindex(dim) ])
vol = vol.reshape(dim)
m = hmesh.volumetric_isocontour(vol, plo, phi, make_triangles=False, high_is_inside=False, dual_connectivity=True)
d = np.clip(D.signed_distance(m.positions()), -0.0005, 0.0005)
v.display(m, mode='w', reset_view=True)
v.display(m, mode='s', data=d)
