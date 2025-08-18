""" DEPRECATED!!!! This module provides a kD-tree implementation but specialized to 3D. The reason why it is deprecated is that scipy.spatial contains a very flexible KDTree class that is much more generic and obviates the need for GEL's KDTree. Going forward, the KDTree may be removed from PyGEL but not from the C++ GEL API. """

from pygel3d import lib_py_gel, Vec3dVector, IntVector
import ctypes as ct

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
    def insert(self,p,v):
        """ Insert v at 3D point given by p. Insert should be called before
        calling build. """
        lib_py_gel.I3DTree_insert(self.obj, p[0],p[1],p[2],v)
    def build(self):
        """ Build the tree. This function call makes the tree searchable. It is
        assumed that all calls to insert come before calling this function."""
        lib_py_gel.I3DTree_build(self.obj)
    def closest_point(self, p, r):
        """ Search for point closest to p within a max radius r. """
        key, val = lib_py_gel.I3DTree_closest_point(self.obj, p[0], p[1], p[2], r)
        if key and len(key) == 3:
            return (list(key), val)
        return None
    def in_sphere(self, p, r):
        """ Retrieve all points within a radius r of p. """
        keys, vals = lib_py_gel.I3DTree_in_sphere(self.obj, p[0], p[1], p[2], r)
        return (list(keys), list(vals))
    def m_closest_points(self, p, r, m):
        """ Retrieve m closest points to p within radius r. """
        keys, vals = lib_py_gel.I3DTree_m_closest_points(self.obj, p[0], p[1], p[2], r, m)
        return (list(keys), list(vals))
