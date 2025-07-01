from numpy import array, zeros, clip, ndindex, ceil, median, cos, sin, meshgrid, linspace
from pygel3d import hmesh, gl_display as gl

v = gl.Viewer()

# First, create an empty volume and polygonize the iso-contour
# (simply checks robustness)
dim = (40,50,60)
vol = zeros(dim)
m = hmesh.volumetric_isocontour(vol)
print("Empty volume has", len(m.vertices()), "vertices and", len(m.faces()), "faces")
v.display(m)

# Second test: create a volume containing a gyroid and display its isocontour
X,Y,Z = meshgrid(linspace(0,dim[0]*0.25, dim[0]),
                 linspace(0,dim[1]*0.25, dim[1]),
                 linspace(0,dim[2]*0.25, dim[2]))
vol = sin(X)*cos(Y)+sin(Y)*cos(Z)+sin(Z)*cos(X)
m = hmesh.volumetric_isocontour(vol)
v.display(m, mode='g', smooth=False, reset_view=True)

# Third test: create a volume containing a gyroid but cut out some slices by setting
# the voxels to NaN and display the isocontour. Surfaces are now open.
vol[:,20:30,:] = float('nan')
vol[15:25,:,:] = float('nan')
vol[:,:,25:35] = float('nan')
m = hmesh.volumetric_isocontour(vol)
v.display(m, mode='g')

# Fourth test: create a volume containing a half space delimited by the volume and 
# display the isocontour. This time the surface is closed.
dim = (16,16,16)
vol = zeros(dim)
n = array((.5,.3,.1))
p0 = array((6,6,6))
for idx in ndindex(dim):
    p = array(idx)
    vol[idx] = n @ (p - p0)
m = hmesh.volumetric_isocontour(vol, make_triangles=False, high_is_inside=False)
v.display(m,smooth=False,reset_view=True)

# Fifth test, create a volume from a signed distance field of the bunny and 
# display the mesh contoured using dual contouring.
m_in = hmesh.load("../../../data/bunny.obj")
D = hmesh.MeshDistance(m_in)
plo, phi = hmesh.bbox(m_in)
dim = array([int(x) for x in ceil((phi-plo) / 0.0025)], dtype=int)
BSZ = 4 # Block size for the distance field 
dim = BSZ * (dim // BSZ)  # Make sure the dimensions are multiples of BSZ
spacing = (phi-plo)/(dim - 1)
vol = D.signed_distance([array(idx)*spacing + plo for idx in ndindex(*dim)])
vol = vol.reshape(dim)
m = hmesh.volumetric_isocontour(vol, plo, phi, high_is_inside=False, dual_connectivity=True)
d = D.signed_distance(m.positions())
v.display(m, mode='s', data=d, reset_view=True)

# Sixth test: compare a simplified mesh to an iso-contour of a downsampled volume.
# Simplify the mesh to 5% of its original number of vertices
hmesh.quadric_simplify(m, 0.05)
v.display(m, mode='w', smooth=False)
# Next, downfilter the volume by averaging BSZ^3 blocks of voxels and polygonize the 
# iso-contour. Then display the mesh with the distance field as a scalar field.
vol_little = vol.reshape((vol.shape[0]//BSZ, BSZ, vol.shape[1]//BSZ, BSZ, vol.shape[2]//BSZ, BSZ))
# Below we downfilter using median, but you can use min, max, mean, etc.
vol_little = median(vol_little, axis=(1,3,5))
m = hmesh.volumetric_isocontour(vol_little, plo, phi, high_is_inside=False, dual_connectivity=True)
v.display(m, mode='w', smooth=False)

