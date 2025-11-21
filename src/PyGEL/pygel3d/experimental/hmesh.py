import numpy as np
import ctypes as ct
from numpy.typing import ArrayLike
from pygel3d.hmesh import Manifold
from pygel3d import lib_py_gel

def rsr_recon(verts: ArrayLike,
              normals: ArrayLike=None,
              use_Euclidean_distance: bool=False,
              genus: int=-1,
              k: int=70,
              r: float=20,
              theta: float=60,
              n: int=50) -> Manifold:
    """ RsR Reconstruction. The first argument, verts, is the point cloud. The next argument,
        normals, are the normals associated with the vertices or empty list (default) if normals
        need to be estimated during reconstruction. use_Euclidean_distance should be true if we
        can use the Euclidean rather than projected distance. Set to true only for noise free
        point clouds. genus is used to constrain the genus of the reconstructed object. genus
        defaults to -1, meaning unknown genus. k is the number of nearest neighbors for each point,
        r is the maximum distance to farthest neighbor measured in multiples of average distance,
        theta is the threshold on angles between normals: two points are only connected if the angle
        between their normals is less than theta. Finally, n is the threshold on the distance between
        vertices that are connected by handle edges (check paper). For large n, it is harder for
        the algorithm to add handles. """
    m = Manifold()
    verts_data = np.asarray(verts, dtype=ct.c_double, order='F')
    n_verts = len(verts)
    n_normal = 0 if normals is None else len(normals)
    if(n_normal==0):
        normals = [[]]
    normal_data = np.asarray(normals, dtype=ct.c_double, order='F')

    lib_py_gel.rsr_recon_experimental(m.obj, verts_data, normal_data, n_verts, n_normal,
                         use_Euclidean_distance, genus, k, r, theta, n)
    return m

def hrsr_recon(verts: ArrayLike,
              normals: ArrayLike=None,
              collapse_iters = 4,
              use_Euclidean_distance: bool=False,
              genus: int=-1,
              k: int=70,
              r: float=20,
              theta: float=60,
              n: int=50,
              skip_reexpansion = False) -> Manifold:
    m = Manifold()
    verts_data = np.asarray(verts, dtype=ct.c_double, order='F')
    n_verts = len(verts)
    n_normal = 0 if normals is None else len(normals)
    if(n_normal==0):
        normals = [[]]
    normal_data = np.asarray(normals, dtype=ct.c_double, order='F')

    lib_py_gel.hrsr_recon_experimental(m.obj, verts_data, normal_data, n_verts, n_normal,
                                      collapse_iters, use_Euclidean_distance, genus, k, r, theta, n, skip_reexpansion)
    return m