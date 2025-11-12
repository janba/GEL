""" DEPRECATED!!!! This module provides a kD-tree implementation but specialized to 3D. The reason why it is deprecated is that scipy.spatial contains a very flexible KDTree class that is much more generic and obviates the need for GEL's KDTree. Going forward, the KDTree may be removed from PyGEL but not from the C++ GEL API. """

import ctypes as ct
from numpy.typing import ArrayLike
from pygel3d import lib_py_gel, Vec3dVector, IntVector

class I3DTree:
    """ kD tree specialized for 3D keys and integer values.
    This tree data structure is useful for storing 3D points and
    associated integer values - typically indices. There is also
    a more general kd tree in scipy.spatial if this one does not
    suit your needs. """
    def __init__(self):
        self.obj = lib_py_gel.I3DTree_new()
    def __del__(self):
        lib_py_gel.I3DTree_delete(self.obj)
    def insert(self,p: ArrayLike, v: int):
        """ Insert v at 3D point given by p. Insert should be called before
        calling build. """
        lib_py_gel.I3DTree_insert(self.obj, p[0],p[1],p[2],v)
    def build(self):
        """ Build the tree. This function call makes the tree searchable. It is
        assumed that all calls to insert come before calling this function."""
        lib_py_gel.I3DTree_build(self.obj)
    def closest_point(self, p: ArrayLike, r: float) -> tuple[list[float], int]|None:
        """ Search for point closest to p within a max radius r.
        This function should only be called after build. """
        key = (ct.c_double * 3)()
        val = ct.c_size_t()
        n = lib_py_gel.I3DTree_closest_point(self.obj, p[0],p[1],p[2],r,ct.byref(key),ct.byref(val))
        if n==1:
            return ([key[0],key[1],key[2]],val.value)
        return None
    def in_sphere(self, p: ArrayLike, r: float) -> tuple[Vec3dVector, IntVector]:
        """ Retrieve all points within a radius r of p.
        This function should only be called after build. """
        keys = Vec3dVector()
        vals = IntVector()
        n = lib_py_gel.I3DTree_in_sphere(self.obj, p[0],p[1],p[2],r,keys.obj,vals.obj)
        return (keys,vals)
