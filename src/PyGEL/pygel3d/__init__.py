""" PyGEL is a collection of classes and functions that expose features in the
GEL library. The primary purpose of PyGEL (and GEL) is to be useful for geometry
processing tasks. Especially tasks that involve 3D polygonal meshes, but there is
also a graph component useful e.g. for skeletonization. The PyGEL package is called
pygel3d and it contains five modules:

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
"""
__all__ = ["hmesh", "graph", "gl_display", "jupyter_display", "spatial"]

import ctypes as ct
import numpy as np
import os
from sys import platform,prefix

def _get_script_path():
    return os.path.dirname(__file__)

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

# The following many lines explitize the arguments and return types of the C API

# IntVector
lib_py_gel.IntVector_new.restype = ct.c_void_p
lib_py_gel.IntVector_get.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.IntVector_size.argtypes = (ct.c_void_p,)
lib_py_gel.IntVector_size.restype = ct.c_size_t
lib_py_gel.IntVector_delete.argtypes = (ct.c_void_p,)


# Vec3dVector
lib_py_gel.Vec3dVector_new.restype = ct.c_void_p
lib_py_gel.Vec3dVector_get.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Vec3dVector_get.restype = ct.POINTER(ct.c_double)
lib_py_gel.Vec3dVector_size.argtypes = (ct.c_void_p,)
lib_py_gel.Vec3dVector_size.restype = ct.c_size_t
lib_py_gel.Vec3dVector_delete.argtypes = (ct.c_void_p,)

# I3DTree
lib_py_gel.I3DTree_new.restype = ct.c_void_p
lib_py_gel.I3DTree_delete.argtypes = (ct.c_void_p,)
lib_py_gel.I3DTree_insert.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_size_t)
lib_py_gel.I3DTree_build.argtypes = (ct.c_void_p,)
lib_py_gel.I3DTree_closest_point.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.POINTER(ct.c_double*3), ct.POINTER(ct.c_size_t))
lib_py_gel.I3DTree_in_sphere.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.c_void_p,ct.c_void_p)

# Manifold class
lib_py_gel.Manifold_from_triangles.argtypes = (ct.c_size_t,ct.c_size_t, np.ctypeslib.ndpointer(ct.c_double), np.ctypeslib.ndpointer(ct.c_int))
lib_py_gel.Manifold_from_triangles.restype = ct.c_void_p
lib_py_gel.Manifold_from_points.argtypes = (ct.c_size_t,np.ctypeslib.ndpointer(ct.c_double), np.ctypeslib.ndpointer(ct.c_double),np.ctypeslib.ndpointer(ct.c_double))
lib_py_gel.Manifold_from_points.restype = ct.c_void_p
lib_py_gel.Manifold_new.restype = ct.c_void_p
lib_py_gel.Manifold_copy.restype = ct.c_void_p
lib_py_gel.Manifold_copy.argtypes = (ct.c_void_p,)
lib_py_gel.Manifold_delete.argtypes = (ct.c_void_p,)
lib_py_gel.Manifold_positions.restype = ct.c_size_t
lib_py_gel.Manifold_positions.argtypes = (ct.c_void_p, ct.POINTER(ct.POINTER(ct.c_double)))
lib_py_gel.Manifold_no_allocated_vertices.restype = ct.c_size_t
lib_py_gel.Manifold_no_allocated_vertices.argtypes = (ct.c_void_p,)
lib_py_gel.Manifold_no_allocated_faces.restype = ct.c_size_t
lib_py_gel.Manifold_no_allocated_faces.argtypes = (ct.c_void_p,)
lib_py_gel.Manifold_no_allocated_halfedges.restype = ct.c_size_t
lib_py_gel.Manifold_no_allocated_halfedges.argtypes = (ct.c_void_p,)
lib_py_gel.Manifold_vertices.restype = ct.c_size_t
lib_py_gel.Manifold_vertices.argtypes = (ct.c_void_p, ct.c_void_p)
lib_py_gel.Manifold_faces.restype = ct.c_size_t
lib_py_gel.Manifold_faces.argtypes = (ct.c_void_p, ct.c_void_p)
lib_py_gel.Manifold_halfedges.restype = ct.c_size_t
lib_py_gel.Manifold_halfedges.argtypes = (ct.c_void_p,ct.c_void_p)
lib_py_gel.Manifold_circulate_vertex.restype = ct.c_size_t
lib_py_gel.Manifold_circulate_vertex.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_char, ct.c_void_p)
lib_py_gel.Manifold_circulate_face.restype = ct.c_size_t
lib_py_gel.Manifold_circulate_face.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_char, ct.c_void_p)
lib_py_gel.Manifold_add_face.argtypes = (ct.c_void_p, ct.c_size_t, np.ctypeslib.ndpointer(ct.c_double))
lib_py_gel.Manifold_remove_face.restype = ct.c_bool
lib_py_gel.Manifold_remove_face.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_remove_edge.restype = ct.c_bool
lib_py_gel.Manifold_remove_edge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_remove_vertex.restype = ct.c_bool
lib_py_gel.Manifold_remove_vertex.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_vertex_in_use.restype = ct.c_bool
lib_py_gel.Manifold_vertex_in_use.argtypes = (ct.c_void_p,ct.c_size_t)
lib_py_gel.Manifold_face_in_use.restype = ct.c_bool
lib_py_gel.Manifold_face_in_use.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_halfedge_in_use.restype = ct.c_bool
lib_py_gel.Manifold_halfedge_in_use.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_flip_edge.restype =  ct.c_bool
lib_py_gel.Manifold_flip_edge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_collapse_edge.restype = ct.c_bool
lib_py_gel.Manifold_collapse_edge.argtypes = (ct.c_void_p,ct.c_size_t,ct.c_bool)
lib_py_gel.Manifold_split_face_by_edge.restype = ct.c_size_t
lib_py_gel.Manifold_split_face_by_edge.argtypes = (ct.c_void_p, ct.c_size_t,ct.c_size_t,ct.c_size_t)
lib_py_gel.Manifold_split_face_by_vertex.restype = ct.c_size_t
lib_py_gel.Manifold_split_face_by_vertex.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_split_edge.restype = ct.c_size_t
lib_py_gel.Manifold_split_edge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Manifold_stitch_boundary_edges.restype = ct.c_bool
lib_py_gel.Manifold_stitch_boundary_edges.argtypes = (ct.c_void_p, ct.c_size_t,ct.c_size_t)
lib_py_gel.Manifold_merge_faces.restype = ct.c_bool
lib_py_gel.Manifold_merge_faces.argtypes = (ct.c_void_p, ct.c_size_t,ct.c_size_t)
lib_py_gel.Manifold_close_hole.argtypes = (ct.c_void_p,ct.c_size_t)
lib_py_gel.Manifold_cleanup.argtypes = (ct.c_void_p,)

# Walker is a helper class assisting us in navigating a mesh.
# Not directly expose in PyGEL3D
lib_py_gel.Walker_next_halfedge.restype = ct.c_size_t
lib_py_gel.Walker_next_halfedge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Walker_prev_halfedge.restype = ct.c_size_t
lib_py_gel.Walker_prev_halfedge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Walker_opposite_halfedge.restype = ct.c_size_t
lib_py_gel.Walker_opposite_halfedge.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Walker_incident_face.restype = ct.c_size_t
lib_py_gel.Walker_incident_face.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Walker_incident_vertex.restype = ct.c_size_t
lib_py_gel.Walker_incident_vertex.argtypes = (ct.c_void_p, ct.c_size_t)

# A list of helper functions
lib_py_gel.is_halfedge_at_boundary.restype = ct.c_bool
lib_py_gel.is_halfedge_at_boundary.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.is_vertex_at_boundary.restype = ct.c_bool
lib_py_gel.is_vertex_at_boundary.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.length.restype = ct.c_double
lib_py_gel.length.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.boundary_edge.restype = ct.c_bool
lib_py_gel.boundary_edge.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_size_t)
lib_py_gel.valency.restype = ct.c_size_t
lib_py_gel.valency.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.vertex_normal.argtypes = (ct.c_void_p, ct.c_size_t, ct.POINTER(ct.c_double*3))
lib_py_gel.connected.restype = ct.c_bool
lib_py_gel.connected.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_size_t)
lib_py_gel.no_edges.restype = ct.c_size_t
lib_py_gel.no_edges.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.face_normal.argtypes = (ct.c_void_p, ct.c_size_t, ct.POINTER(ct.c_double*3))
lib_py_gel.area.restype = ct.c_double
lib_py_gel.area.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.perimeter.restype = ct.c_double
lib_py_gel.perimeter.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.centre.argtypes = (ct.c_void_p, ct.c_size_t, ct.POINTER(ct.c_double*3))
lib_py_gel.valid.restype = ct.c_bool
lib_py_gel.valid.argtypes = (ct.c_void_p,)
lib_py_gel.closed.restype = ct.c_bool
lib_py_gel.closed.argtypes = (ct.c_void_p,)
lib_py_gel.bbox.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double*3),ct.POINTER(ct.c_double*3))
lib_py_gel.bsphere.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double*3), ct.POINTER(ct.c_double))
lib_py_gel.stitch_mesh.argtypes = (ct.c_void_p,ct.c_double)
lib_py_gel.stitch_mesh.restype = ct.c_int
lib_py_gel.obj_save.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.off_save.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.x3d_save.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.obj_load.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.off_load.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.ply_load.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.x3d_load.argtypes = (ct.c_char_p, ct.c_void_p)
lib_py_gel.remove_caps.argtypes = (ct.c_void_p, ct.c_float)
lib_py_gel.remove_needles.argtypes = (ct.c_void_p, ct.c_float, ct.c_bool)
lib_py_gel.close_holes.argtypes = (ct.c_void_p,ct.c_int)
lib_py_gel.flip_orientation.argtypes = (ct.c_void_p,)
lib_py_gel.merge_coincident_boundary_vertices.argtypes = (ct.c_void_p, ct.c_double)
lib_py_gel.minimize_curvature.argtypes = (ct.c_void_p,ct.c_bool)
lib_py_gel.minimize_dihedral_angle.argtypes = (ct.c_void_p, ct.c_int, ct.c_bool, ct.c_bool, ct.c_double)
lib_py_gel.maximize_min_angle.argtypes = (ct.c_void_p,ct.c_float,ct.c_bool)
lib_py_gel.optimize_valency.argtypes = (ct.c_void_p,ct.c_bool)
lib_py_gel.randomize_mesh.argtypes = (ct.c_void_p,ct.c_int)
lib_py_gel.quadric_simplify.argtypes = (ct.c_void_p,ct.c_double,ct.c_double,ct.c_bool)
lib_py_gel.average_edge_length.argtypes = (ct.c_void_p,)
lib_py_gel.average_edge_length.restype = ct.c_float
lib_py_gel.median_edge_length.argtypes = (ct.c_void_p,)
lib_py_gel.median_edge_length.restype = ct.c_float
lib_py_gel.refine_edges.argtypes = (ct.c_void_p,ct.c_float)
lib_py_gel.refine_edges.restype = ct.c_int
lib_py_gel.cc_split.argtypes = (ct.c_void_p,)
lib_py_gel.loop_split.argtypes = (ct.c_void_p,)
lib_py_gel.root3_subdivide.argtypes = (ct.c_void_p,)
lib_py_gel.rootCC_subdivide.argtypes = (ct.c_void_p,)
lib_py_gel.butterfly_subdivide.argtypes = (ct.c_void_p,)
lib_py_gel.cc_smooth.argtypes = (ct.c_void_p,)
lib_py_gel.loop_smooth.argtypes = (ct.c_void_p,)
lib_py_gel.ear_clip_triangulate.argtypes = (ct.c_void_p,)
lib_py_gel.shortest_edge_triangulate.argtypes = (ct.c_void_p,)
lib_py_gel.graph_to_feq.argtypes = (ct.c_void_p, ct.c_void_p, ct.POINTER(ct.c_double))
lib_py_gel.graph_to_feq.restype = ct.c_void_p


# MeshDistance allows us to compute the signed distance to a mesh
lib_py_gel.MeshDistance_new.restype = ct.c_void_p
lib_py_gel.MeshDistance_new.argtypes = (ct.c_void_p,)
lib_py_gel.MeshDistance_signed_distance.argtypes = (ct.c_void_p,ct.c_int, ct.POINTER(ct.c_float),ct.POINTER(ct.c_float),ct.c_float)
lib_py_gel.MeshDistance_ray_inside_test.argtypes = (ct.c_void_p,ct.c_int, ct.POINTER(ct.c_float),ct.POINTER(ct.c_int),ct.c_int)
lib_py_gel.MeshDistance_delete.argtypes = (ct.c_void_p,)

# The Graph class
lib_py_gel.Graph_new.restype = ct.c_void_p
lib_py_gel.Graph_copy.restype = ct.c_void_p
lib_py_gel.Graph_copy.argtypes = (ct.c_void_p,)
lib_py_gel.Graph_delete.argtypes = (ct.c_void_p,)
lib_py_gel.Graph_clear.argtypes = (ct.c_void_p,)
lib_py_gel.Graph_cleanup.argtypes = (ct.c_void_p,)
lib_py_gel.Graph_nodes.argtypes = (ct.c_void_p, ct.c_void_p)
lib_py_gel.Graph_nodes.restype = ct.c_size_t
lib_py_gel.Graph_neighbors.restype = ct.c_size_t
lib_py_gel.Graph_neighbors.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_void_p, ct.c_char)
lib_py_gel.Graph_positions.argtypes = (ct.c_void_p,ct.POINTER(ct.POINTER(ct.c_double)))
lib_py_gel.Graph_positions.restype = ct.c_size_t
lib_py_gel.Graph_average_edge_length.argtypes = (ct.c_void_p,)
lib_py_gel.Graph_add_node.argtypes = (ct.c_void_p, np.ctypeslib.ndpointer(ct.c_double))
lib_py_gel.Graph_add_node.restype = ct.c_size_t
lib_py_gel.Graph_remove_node.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Graph_node_in_use.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Graph_connect_nodes.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_size_t)
lib_py_gel.Graph_connect_nodes.restype = ct.c_size_t
lib_py_gel.Graph_disconnect_nodes.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_size_t)
lib_py_gel.Graph_merge_nodes.argtypes = (ct.c_void_p, ct.c_size_t, ct.c_size_t, ct.c_bool)

# Graph functions
lib_py_gel.graph_from_mesh.argtypes = (ct.c_void_p, ct.c_void_p)
lib_py_gel.graph_load.argtypes = (ct.c_void_p, ct.c_char_p)
lib_py_gel.graph_load.restype = ct.c_void_p
lib_py_gel.graph_save.argtypes = (ct.c_void_p, ct.c_char_p)
lib_py_gel.graph_save.restype = ct.c_bool
lib_py_gel.graph_to_mesh_cyl.argtypes = (ct.c_void_p, ct.c_void_p, ct.c_float)
lib_py_gel.graph_to_mesh_cyl.restype = ct.c_void_p
lib_py_gel.graph_smooth.argtypes = (ct.c_void_p, ct.c_int, ct.c_float)
lib_py_gel.graph_edge_contract.argtypes = (ct.c_void_p, ct.c_double)
lib_py_gel.graph_prune.argtypes = (ct.c_void_p,)
lib_py_gel.graph_LS_skeleton.argtypes = (ct.c_void_p, ct.c_void_p, ct.c_void_p, ct.c_bool)
lib_py_gel.graph_front_skeleton.argtypes = (ct.c_void_p, ct.c_void_p, ct.c_void_p, ct.c_int, ct.POINTER(ct.c_double))

class IntVector:
    """ Vector of integer values.
    This is a simple class that implements iteration and index based
    retrieval. Allocation happens in a call to libPyGEL. Since memory
    is managed by the PyGEL library, the vector can be resized by library
    functions. Generally not used directly by PyGEL3D users."""
    def __init__(self):
        self.obj = lib_py_gel.IntVector_new(0)
    def __del__(self):
        lib_py_gel.IntVector_delete(self.obj)
    def __len__(self):
        return int(lib_py_gel.IntVector_size(self.obj))
    def __getitem__(self,key):
        return lib_py_gel.IntVector_get(self.obj,key)
    def __iter__(self):
        n = lib_py_gel.IntVector_size(self.obj)
        for i in range(0,n):
            yield lib_py_gel.IntVector_get(self.obj, i)

class Vec3dVector:
    """ Vector of 3D vectors.
    This is a simple class that implements iteration and index based
    retrieval. Allocation happens in a call to libPyGEL. Since memory
    is managed by the PyGEL library, the vector can be resized by library
    functions. nerally not used directly by PyGEL3D users."""
    def __init__(self):
        self.obj = lib_py_gel.Vec3dVector_new(0)
    def __del__(self):
        lib_py_gel.Vec3dVector_delete(self.obj)
    def __len__(self):
        return int(lib_py_gel.Vec3dVector_size(self.obj))
    def __getitem__(self,key):
        return lib_py_gel.Vec3dVector_get(self.obj,key)
    def __iter__(self):
        n = lib_py_gel.Vec3dVector_size(self.obj)
        for i in range(0,n):
            data = lib_py_gel.Vec3dVector_get(self.obj, i)
            yield [data[0], data[1], data[2]]
