from pygel3d import hmesh
from sys import argv
from numpy import array, median
import polyscope as ps

ps.init()
ps.set_ground_plane_mode("shadow_only")

def normalize(v):
    norm = (v @ v) ** 0.5
    if norm > 0:
        return v / norm
    return v

def smooth_vector_field(m: hmesh.Manifold, data, niter=10):
    for _ in range(niter):
        data_new = data.copy()
        for v in m.vertices():
            nbors =  m.circulate_vertex(v)
            nsign = [ 1 if data[n]@data[v] > 0 else -1 for n in nbors ]
            ndata = [ data[n]*nsign[i] for i,n in enumerate(nbors) ]
            avg = sum(ndata) / len(ndata)
            data_new[v] = normalize(0.5 * data[v] + 0.5 * avg)
        data = data_new
    return data

fn = "../../../data/Solids/torus.obj" if len(argv) < 2 else argv[1]
m = hmesh.load(fn)
hmesh.triangulate(m)

verts = m.positions()
faces = [ list(m.circulate_face(f)) for f in m.faces() ]
ps.register_surface_mesh(fn, verts, faces)

pc = [ m.principal_curvatures(v) for v in m.vertices() ]
kmin = array([ pc[v][0] for v in m.vertices()])
ps.get_surface_mesh(fn).add_scalar_quantity("min curvature", kmin, defined_on='vertices', enabled=True)
dirmin = smooth_vector_field(m, array([ pc[v][2] for v in m.vertices()]), niter=10)
ps.get_surface_mesh(fn).add_vector_quantity("min curvature direction", dirmin, defined_on='vertices', enabled=True)

kmax = array([ pc[v][1] for v in m.vertices()])
ps.get_surface_mesh(fn).add_scalar_quantity("max curvature", kmax, defined_on='vertices', enabled=False)
dirmax = smooth_vector_field(m, array([ pc[v][3] for v in m.vertices()]), niter=10)
ps.get_surface_mesh(fn).add_vector_quantity("max curvature direction", dirmax, defined_on='vertices', enabled=False)

K = array([ m.gaussian_curvature(v) for v in m.vertices()])
K2 = kmin * kmax
ps.get_surface_mesh(fn).add_scalar_quantity("Gaussian curvature", K, defined_on='vertices', enabled=False)
ps.get_surface_mesh(fn).add_scalar_quantity("Gaussian curvature (product of principal)", K2, defined_on='vertices', enabled=False)

H = array([ m.mean_curvature(v) for v in m.vertices() ])
H2 = (kmin + kmax) / 2
ps.get_surface_mesh(fn).add_scalar_quantity("Mean curvature", H, defined_on='vertices', enabled=False)
ps.get_surface_mesh(fn).add_scalar_quantity("Mean curvature (average of principal)", H2, defined_on='vertices', enabled=False)

ps.show()
