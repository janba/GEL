from pygel3d import hmesh, gl_display as gl
from sys import argv
from numpy import array
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fn = "../../../data/ReferenceMeshes/torus.obj" if len(argv) < 2 else argv[1]
m = hmesh.load(fn)

v = gl.Viewer()

pc = [ m.principal_curvatures(v) for v in m.vertices() ]
kmin = array([ pc[v][0] for v in m.vertices()])
kmax = array([ pc[v][1] for v in m.vertices()])
print("min curvature range: ", kmin.min(), kmin.max())
print("max curvature range: ", kmax.min(), kmax.max())
v.display(m, mode='s', data=kmin)
v.display(m, mode='s', data=kmax)
dirmin = array([ pc[v][2] for v in m.vertices()])
dirmax = array([ pc[v][3] for v in m.vertices()])

K2 = kmin * kmax
H2 = (kmin + kmax) / 2

K = array([ m.gaussian_curvature(v) for v in m.vertices()])
v.display(m, mode='s', data=K)
v.display(m, mode='s', data=K2)
print("Gaussian curvature range: ", K.min(), K.max())
print("Gaussian curvature range (2): ", K2.min(), K2.max())

H = array([ m.mean_curvature(v) for v in m.vertices() ])
v.display(m, mode='s', data=H)
v.display(m, mode='s', data=H2)
print("Mean curvature range: ", H.min(), H.max())
print("Mean curvature range (2): ", H2.min(), H2.max())

# Plot dirmin vectors
pos = m.positions()
positions = array([pos[v] for v in m.vertices()])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(positions[:,0], positions[:,1], positions[:,2],
          dirmin[:,0], dirmin[:,1], dirmin[:,2], length=0.1, normalize=True, color='b')
ax.set_title('dirmin vectors')
plt.show()

# Plot dirmax vectors
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(positions[:,0], positions[:,1], positions[:,2],
          dirmax[:,0], dirmax[:,1], dirmax[:,2], length=0.1   , normalize=True, color='r')
ax.set_title('dirmax vectors')
plt.show()