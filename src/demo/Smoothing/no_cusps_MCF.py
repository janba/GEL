from itertools import permutations
from sys import argv
from numpy import zeros, cross, average
from numpy.linalg import norm, det
from scipy.linalg import solve, cho_factor, cho_solve
from scipy.sparse.linalg import cg
from pygel3d import hmesh, gl_display as gl
#from cvxopt import matrix, spmatrix, cholmod

def normalize(v):
    return v / norm(v)

def vertex_barycenter(m):
    pos = m.positions()
    return average(pos, axis=0)

def rescale(m, vol_target, c_target):
    vol_new = hmesh.volume(m)
    s = (vol_target/vol_new)**(1.0/3.0)
    c = vertex_barycenter(m)
    for i in m.vertices():
        pos[i] = s*(pos[i]-c) + c_target
    return vol_new


# Load the mesh and make safe copy of vertex positions
# m = hmesh.load("arma.obj")
m = hmesh.load(argv[1])
pos = m.positions()
N = m.no_allocated_vertices()
v = gl.Viewer()

v.display(m,smooth=False, once=False) # HIT ESC IN THE GRAPHICS WINDOW TO CONTINUE !!!!

print("Making cotan Laplacian")
ael = hmesh.average_edge_length(m)
L = create_stiffness_matrix(m)
print("Done with cotan Laplacian")
vol0 = hmesh.volume(m)
c0 = vertex_barycenter(m)
N = 250
for iter in range(N):
    M = create_mass_matrix(m)
    dt = ael * 10 * ((iter+1)/N)**0.5
    # Solve using Cholesky
    pos[:] = cho_solve(cho_factor(M-dt*L), M@pos, check_finite=False)

    # Solve using conventional solver
    # pos[:] = solve(M-L, M@pos)

    # Solve using CG
    # B = M@pos
    # K = M-L
    # pos[:,0],_ = cg(K, B[:,0])
    # pos[:,1],_ = cg(K, B[:,1])
    # pos[:,2],_ = cg(K, B[:,2])
    rescale(m, vol0, c0)
    v.display(m,smooth=False, once=True)

new_vol = hmesh.volume(m)
print("Volume change: ", new_vol/vol0)

del v
