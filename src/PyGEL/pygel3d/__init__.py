""" PyGEL is a collection of classes and functions for geometry processing tasks. 
Especially tasks that involve 3D polygonal meshes, but there is also a graph component 
useful e.g. for skeletonization. The PyGEL package is called pygel3d and it contains 
five modules:

hmesh provides Manifold which is a class that represents polygonal meshes using the
halfedge representation. hmesh also provides a slew of functions for manipulating
polygonal meshes and the MeshDistance class which makes it simple to compute the
distance to a triangle mesh.

graph contains the Graph class which is used for graphs: i.e. collections of
vertices (in 3D) connected by edges. Unlike a Manifold, a Graph does not have to
represent a surface. There are also some associated functions which may be useful:
in particular, there is the LS_skeletonize function which computes a curve skeleton
from a Graph and returns the result as a new Graph.

gl_display provides the Viewer class which makes it simple to visualize meshes and
graphs.

jupyter_display makes it easy to use PyGEL in the context of a Jupyter Notebook.
This module contains a function that allows you to create a widget for interactively
visualizing a mesh or a graph in a Notebook. The feature is based on the Plotly
library and it is possible to export the resulting notebooks to HTML while preserving
the interactive 3D graphics in the notebook.

spatial contains the I3DTree class which is simply a kD-tree specialized for mapping
3D points to integers - typically indices. Of course, scipy.spatial has a more
generic class, so this is perhaps not the most important part of PyGEL.

PyGEL is based on the C++ GEL library and provides a Python interface for most but not
all of the functionality of GEL. 
"""
__all__ = ["hmesh", "graph", "gl_display", "jupyter_display", "spatial"]

import os
import numpy as np

def _get_script_path():
    return os.path.dirname(__file__)

try:
    # Try to import the pybind11 PyGEL module
    from . import PyGEL as lib_py_gel
    # Set InvalidIndex from the module
    InvalidIndex = lib_py_gel.InvalidIndex
except ImportError:
    # Fallback to ctypes if pybind11 module not available
    import ctypes as ct
    from numpy.ctypeslib import ndpointer
    from sys import platform
    
    def _get_lib_name():
        if platform == "darwin":
            return "libPyGEL.dylib"
        if platform == "win32":
            return "PyGEL.dll"
        return "libPyGEL.so"

    # Load PyGEL the Python GEL bridge library
    lib_py_gel = ct.cdll.LoadLibrary(_get_script_path() + "/" + _get_lib_name())

    # An InvalidIndex is just a special integer value.
    InvalidIndex = ct.c_size_t.in_dll(lib_py_gel, "InvalidIndex").value
    
    # Note: Ctypes definitions would go here, but we're using pybind11 now

class IntVector:
    """ Vector of integer values.
    This is a simple class that implements iteration and index based
    retrieval. For pybind11, this will be a thin wrapper around Python lists."""
    def __init__(self, data=None):
        if hasattr(lib_py_gel, 'IntVector_new'):
            # pybind11 version
            self.obj = lib_py_gel.IntVector_new(0) if data is None else data
        else:
            # fallback
            self.obj = []
            
    def __del__(self):
        if hasattr(lib_py_gel, 'IntVector_delete') and hasattr(self, 'obj'):
            lib_py_gel.IntVector_delete(self.obj)
            
    def __len__(self):
        if hasattr(lib_py_gel, 'IntVector_size'):
            return int(lib_py_gel.IntVector_size(self.obj))
        return len(self.obj)
        
    def __getitem__(self, key):
        if hasattr(lib_py_gel, 'IntVector_get'):
            return lib_py_gel.IntVector_get(self.obj, key)
        return self.obj[key]
        
    def __iter__(self):
        n = len(self)
        for i in range(n):
            yield self[i]

class Vec3dVector:
    """ Vector of 3D vectors.
    This is a simple class that implements iteration and index based
    retrieval. For pybind11, this will be a thin wrapper."""
    def __init__(self, data=None):
        if hasattr(lib_py_gel, 'Vec3dVector_new'):
            # pybind11 version  
            self.obj = lib_py_gel.Vec3dVector_new(0) if data is None else data
        else:
            # fallback
            self.obj = []
            
    def __del__(self):
        if hasattr(lib_py_gel, 'Vec3dVector_delete') and hasattr(self, 'obj'):
            lib_py_gel.Vec3dVector_delete(self.obj)
            
    def __len__(self):
        if hasattr(lib_py_gel, 'Vec3dVector_size'):
            return int(lib_py_gel.Vec3dVector_size(self.obj))
        return len(self.obj)
        
    def __getitem__(self, key):
        if hasattr(lib_py_gel, 'Vec3dVector_get'):
            return lib_py_gel.Vec3dVector_get(self.obj, key)
        return self.obj[key]
        
    def __iter__(self):
        n = len(self)
        for i in range(n):
            yield self[i]
