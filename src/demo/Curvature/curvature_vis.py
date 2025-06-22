from pygel3d import hmesh, gl_display as gl
from sys import argv
from numpy import array
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def smooth_vector_field(m: hmesh.Manifold, data, niter=10):
    """Smooth a vector field on the mesh using Laplacian smoothing."""
    for _ in range(niter):
        data_new = data.copy()
        for v in m.vertices():
            ndata = [ data[n]*(1 if data[n]@data[v] > 0 else -1) for n in m.circulate_vertex(v) ]
            avg = sum(ndata) / len(ndata)
            data_new[v] = 0.5 * data[v] + 0.5 * avg
        data = data_new
    return data

fn = "../../../data/Solids/sphere.obj" if len(argv) < 2 else argv[1]
m = hmesh.load(fn)
hmesh.triangulate(m)

v = gl.Viewer()

pc = [ m.principal_curvatures(v) for v in m.vertices() ]
kmin = array([ pc[v][0] for v in m.vertices()])
kmax = array([ pc[v][1] for v in m.vertices()])
print("min curvature range: ", kmin.min(), kmin.max())
print("max curvature range: ", kmax.min(), kmax.max())

print("Displaying minimum principal curvature")
v.display(m, mode='s', data=kmin)

print("Displaying maximum principal curvature")
v.display(m, mode='s', data=kmax)

dirmin = smooth_vector_field(m, array([ pc[v][2] for v in m.vertices()]))
print("Displaying minimum principal curvature direction")
v.display(m, mode='l', data=dirmin)

print("Displaying maximum principal curvature direction")
dirmax = smooth_vector_field(m, array([ pc[v][3] for v in m.vertices()]))
v.display(m, mode='l', data=dirmax)

K2 = kmin * kmax
H2 = (kmin + kmax) / 2
K = array([ m.gaussian_curvature(v) for v in m.vertices()])

print("Displaying Gaussian curvature")
v.display(m, mode='s', data=K)

print("Displaying Gaussian curvature (product of principal curvatures)")
v.display(m, mode='s', data=K2)
print("Gaussian curvature range: ", K.min(), K.max())
print("Gaussian curvature range (product of principal curvatures): ", K2.min(), K2.max())

H = array([ m.mean_curvature(v) for v in m.vertices() ])

print("Displaying Mean curvature")
v.display(m, mode='s', data=H)

print("Displaying Mean curvature (average of principal curvatures)")
v.display(m, mode='s', data=H2)
print("Mean curvature range: ", H.min(), H.max())
print("Mean curvature range (average of principal curvatures): ", H2.min(), H2.max())
