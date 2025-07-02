from pygel3d import hmesh, gl_display as gl, graph
import numpy as np
from numba import jit, njit
from numpy.linalg import norm
from random import random
array = np.array

@njit
def random_point_in_unit_sphere():
    while True:
        x = array((random(), random(), random()))
        x *= 2 - 1
        if norm(x) < 1.0:
            return x
    return None

@njit
def find_closest(pos: np.ndarray, x):
    idx = np.argmin([norm(pos[i]-x) for i in g.nodes()])
    return idx 

viewer = gl.Viewer()
g = graph.Graph()
g.add_node([.0,.0,.0])

for _ in range(100):
    pos = g.positions()
    x = random_point_in_unit_sphere()
    for inner_iter in range(100000):
        v = 0.01 * random_point_in_unit_sphere()
        x = x + v
        l = norm(x)
        if l > 1.0:
            x_norm = x/l
            x -= 2*x_norm
        i = find_closest(g.positions(), x)
        if norm(x-pos[i]) < 0.01:
            j = g.add_node(x)
            g.connect_nodes(i, j)
            break
viewer.display(g)



