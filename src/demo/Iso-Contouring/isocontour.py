from math import cos, sin
from numpy import array
import numpy as np
from pygel3d import hmesh, gl_display as gl

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
# display the mesh contoured using dual contouring.
m_in = hmesh.load("../../../data/bunny.obj")
D = hmesh.MeshDistance(m_in)
plo, phi = hmesh.bbox(m_in)
dim = array([int(x) for x in np.ceil((phi-plo) / 0.0025)], dtype=int)
BSZ = 4 # Block size for the distance field 
dim += BSZ - dim % BSZ  # Make sure the dimensions are multiples of BSZ
spacing = (phi-plo)/(dim - 1)
vol = D.signed_distance([array(idx)*spacing + plo for idx in np.ndindex(*dim)])
vol = vol.reshape(dim)
m = hmesh.volumetric_isocontour(vol, plo, phi, high_is_inside=False)
# Simplify the mesh to 5% of its original number of vertices
hmesh.quadric_simplify(m, 0.05)
# Improve the mesh quality by minimizing dihedral angles and maximizing minimum angles
hmesh.minimize_dihedral_angle(m) # Make the mesh geometry better (i.e. reduce curvature)
hmesh.maximize_min_angle(m, dihedral_thresh=0.98) # improve the triangle quality (i.e. more equilateral)
v.display(m, mode='w', reset_view=True, smooth=False)

# Fifth test, downfilter the volume by averaging BSZ^3 blocks of voxels and polygonize the 
# iso-contour. Then display the mesh with the distance field as a scalar field.
vol_little = vol.reshape((vol.shape[0]//BSZ, BSZ, vol.shape[1]//BSZ, BSZ, vol.shape[2]//BSZ, BSZ))
# Below we downfilter using median, but you can use min, max, mean, etc.
vol_little = np.median(vol_little, axis=(1,3,5))
m = hmesh.volumetric_isocontour(vol_little, plo, phi, high_is_inside=False, dual_connectivity=True)
# Show the distance field as a scalar field on the mesh.
d = np.clip(D.signed_distance(m.positions()), -0.005, 0.005)
v.display(m, mode='s', data=d)

# Sixth test, we redisplay the mesh and polygonize without dual connectivity and display the mesh
# again. Note that dual connectivity tends to produce better quality triangles.
v.display(m, mode='w', smooth=False)
m = hmesh.volumetric_isocontour(vol_little, plo, phi, high_is_inside=False, dual_connectivity=False)
v.display(m, mode='w', smooth=False)
