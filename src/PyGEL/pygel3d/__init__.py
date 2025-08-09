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
import sys
from glob import glob
import importlib
import importlib.util
import numpy as np  # kept for API compatibility; modules import this symbol

def _get_script_path():
    return os.path.dirname(__file__)

def _import_pybind_from_package_dir():
    try:
        module_name = "PyGEL"
        module_path = glob(_get_script_path() + "/" + module_name + "*.so")[0]
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        PyGEL = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(PyGEL)
        return PyGEL
    except Exception as e:
        print(f"Error importing {module_name}: {e}")
    raise ImportError(
        "Could not locate the PyGEL pybind11 extension next to pygel3d. "
        "Expected a file like 'PyGEL*.so' (Unix/macOS) or 'PyGEL*.pyd' (Windows) "
        "in the same directory as pygel3d/__init__.py."
    )

# Import only the pybind11 extension located alongside this package in site-packages.
lib_py_gel = _import_pybind_from_package_dir()

# Expose InvalidIndex if provided by the extension
try:
    InvalidIndex = lib_py_gel.InvalidIndex
except AttributeError:
    InvalidIndex = None

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
