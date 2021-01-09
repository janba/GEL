""" PyGEL.gel is a collection of classes and functions that expose features in the
GEL library. The primary purpose of PyGEL (and GEL) is to be useful for geometry
processing tasks. Especially tasks that involve 3D polygonal meshes.

Since the representation of meshes is halfedge based, we are not restricted to
triangles. On the other hand, the library must keep the connectivity information
up to date.

The principal features of PyGEL are

1. Manifold class: a halfedge mesh representation.
2. GLManifoldViewer class: an OpenGL based viewer for Manifold.
3. I3DTree class: a kD-tree data structure for 3D point clouds.
4. MeshDistance class: a data structure for distance queries.
5. A fairly rich library of functions are available for the
Manifold class, allowing users to simplify, optimize, triangulate, refine,
or otherwise manipulate general polygon meshes.
"""

import ctypes as ct
import numpy as np
from numpy.linalg import norm
import os
from sys import platform,prefix

def get_script_path():
    return os.path.dirname(__file__)

def get_lib_name():
    if platform == "darwin":
        return "libPyGEL.dylib"
    if platform == "win32":
        return "PyGEL.dll"
    return "libPyGEL.so"

# Load PyGEL the Python GEL bridge library
lib_py_gel = ct.cdll.LoadLibrary(get_script_path() + "/" + get_lib_name())

InvalidIndex = ct.c_size_t.in_dll(lib_py_gel, "InvalidIndex").value

lib_py_gel.IntVector_new.restype = ct.c_void_p
lib_py_gel.IntVector_get.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.IntVector_size.argtypes = (ct.c_void_p,)
lib_py_gel.IntVector_size.restype = ct.c_size_t
lib_py_gel.IntVector_delete.argtypes = (ct.c_void_p,)
class IntVector:
    """ Vector of integer values.
    This is a simple class that implements iteration and index based
    retrieval. Allocation happens in a call to libPyGEL. Since memory
    is managed by the PyGEL library, the vector can be resized by library
    functions."""
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

lib_py_gel.Vec3dVector_new.restype = ct.c_void_p
lib_py_gel.Vec3dVector_get.argtypes = (ct.c_void_p, ct.c_size_t)
lib_py_gel.Vec3dVector_get.restype = ct.POINTER(ct.c_double)
lib_py_gel.Vec3dVector_size.argtypes = (ct.c_void_p,)
lib_py_gel.Vec3dVector_size.restype = ct.c_size_t
lib_py_gel.Vec3dVector_delete.argtypes = (ct.c_void_p,)
class Vec3dVector:
    """ Vector of 3D vectors.
    This is a simple class that implements iteration and index based
    retrieval. Allocation happens in a call to libPyGEL. Since memory
    is managed by the PyGEL library, the vector can be resized by library
    functions."""
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


lib_py_gel.I3DTree_new.restype = ct.c_void_p
lib_py_gel.I3DTree_delete.argtypes = (ct.c_void_p,)
lib_py_gel.I3DTree_insert.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_size_t)
lib_py_gel.I3DTree_build.argtypes = (ct.c_void_p,)
lib_py_gel.I3DTree_closest_point.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.POINTER(ct.c_double*3), ct.POINTER(ct.c_size_t))
lib_py_gel.I3DTree_in_sphere.argtypes = (ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.c_void_p,ct.c_void_p)
class I3DTree:
    """ kD tree specialized for 3D keys and integer values.
    This tree data structure is useful for storing 3D points and
    associated integer values - typically indices. """
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
    def closest_point(self,p,r):
        """ Search for point closest to p within a max radius r.
        This function should only be called after build. """
        key = (ct.c_double * 3)()
        val = ct.c_size_t()
        n = lib_py_gel.I3DTree_closest_point(self.obj, p[0],p[1],p[2],r,ct.byref(key),ct.byref(val))
        if n==1:
            return ([key[0],key[1],key[2]],val.value)
        return None
    def in_sphere(self, p, r):
        """ Retrieve all points within a radius r of p.
        This function should only be called after build. """
        keys = Vec3dVector()
        vals = IntVector()
        n = lib_py_gel.I3DTree_in_sphere(self.obj, p[0],p[1],p[2],r,keys.obj,vals.obj)
        return (keys,vals)

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

class Manifold:
    """ The Manifold class represents a halfedge based mesh.
    Generally, meshes based on the halfedge representation must be manifold. Hence,
    the name. This class contains methods for mesh manipulation and inspection.
    """
    def __init__(self,orig=None):
        if orig == None:
            self.obj = lib_py_gel.Manifold_new()
        else:
            self.obj = lib_py_gel.Manifold_copy(orig.obj)
    @classmethod
    def from_triangles(cls,vertices, faces):
        m = cls()
        m.obj = lib_py_gel.Manifold_from_triangles(len(vertices),len(faces),np.array(vertices,dtype=np.float64), np.array(faces,dtype=ct.c_int))
        return m
    @classmethod
    def from_points(cls,pts,xaxis=np.array([1,0,0]),yaxis=np.array([0,1,0])):
        m = cls()
        m.obj = lib_py_gel.Manifold_from_points(len(pts),np.array(pts,dtype=np.float64), np.array(xaxis,dtype=np.float64), np.array(yaxis,dtype=np.float64))
        return m
    def __del__(self):
        lib_py_gel.Manifold_delete(self.obj)
    def add_face(self,pts):
        """ Add a face to the Manifold.
        This function takes a list of points as argument and creates a face
        in the mesh with those points as vertices.
        """
        lib_py_gel.Manifold_add_face(self.obj, len(pts), np.array(pts))
    def positions(self):
        """ Retrieve an array containing vertex positions. """
        pos = ct.POINTER(ct.c_double)()
        n = lib_py_gel.Manifold_positions(self.obj, ct.byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
    def no_allocated_vertices(self):
        """ Number of vertices.
        This number could be higher than the number of actually
        used vertices, but corresponds to the number of slots
        allocated."""
        return lib_py_gel.Manifold_no_allocated_vertices(self.obj)
    def no_allocated_faces(self):
        """ Number of faces.
        This number could be higher than the number of actually
        used faces, but corresponds to the number of slots
        allocated."""
        return lib_py_gel.Manifold_no_allocated_faces(self.obj)
    def no_allocated_halfedges(self):
        """ Number of halfedges.
        This number could be higher than the number of actually
        used halfedges, but corresponds to the number of slots
        allocated."""
        return lib_py_gel.Manifold_no_allocated_halfedges(self.obj)
    def vertices(self):
        """ Returns an iterable containing all vertex indices"""
        verts = IntVector()
        n = lib_py_gel.Manifold_vertices(self.obj, verts.obj)
        return verts
    def faces(self):
        """ Returns an iterable containing all face indices"""
        faces = IntVector()
        n = lib_py_gel.Manifold_faces(self.obj, faces.obj)
        return faces
    def halfedges(self):
        """ Returns an iterable containing all halfedge indices"""
        hedges = IntVector()
        n = lib_py_gel.Manifold_halfedges(self.obj, hedges.obj)
        return hedges
    def circulate_vertex(self,v,mode='v'):
        """ Circulate a vertex. Passed a vertex index, v, and second argument,
        mode='f', this function will return an iterable with all faces incident
        on v arranged in counter clockwise order. Similarly, if mode is 'h',
        incident halfedges (outgoing) are returned, and for mode = 'v', all
        neighboring vertices are returned. """
        nbrs = IntVector()
        n = lib_py_gel.Manifold_circulate_vertex(self.obj, v, ct.c_char(mode.encode('ascii')), nbrs.obj)
        return nbrs
    def circulate_face(self,f,mode='v'):
        """ Circulate a face. Passed a face index, f, and second argument,
        mode='f', this function will return an iterable with all faces that
        share an edge with f (in counter clockwise order). If the argument is
        mode='h', the halfedges themselves are returned. For mode='v', the
        incident vertices of the face are returned. """
        nbrs = IntVector()
        n = lib_py_gel.Manifold_circulate_face(self.obj, f, ct.c_char(mode.encode('ascii')), nbrs.obj)
        return nbrs
    def next_halfedge(self,hid):
        """ Returns next halfedge to the one passed as argument. """
        return lib_py_gel.Walker_next_halfedge(self.obj, hid)
    def prev_halfedge(self,hid):
        """ Returns previous halfedge to the one passed as argument. """
        return lib_py_gel.Walker_prev_halfedge(self.obj, hid)
    def opposite_halfedge(self,hid):
        """ Returns opposite halfedge to the one passed as argument. """
        return lib_py_gel.Walker_opposite_halfedge(self.obj, hid)
    def incident_face(self,hid):
        """ Returns face corresponding to halfedge passed as argument. """
        return lib_py_gel.Walker_incident_face(self.obj, hid)
    def incident_vertex(self,hid):
        """ Returns vertex corresponding to halfedge passed as argument. """
        return lib_py_gel.Walker_incident_vertex(self.obj, hid)
    def remove_vertex(self,vid):
        """ Remove a vertex from the Manifold. This function merges all faces
        around the vertex into one and then removes this resulting face. """
        return lib_py_gel.Manifold_remove_vertex(self.obj, vid)
    def remove_face(self,fid):
        """ Removes a face from the Manifold. If it is an interior face it is
        simply replaced by an invalid index. If the face contains boundary
        edges, these are removed. Situations may arise where the mesh is no
        longer manifold because the situation at a boundary vertex is not
        homeomorphic to a half disk. This, we can probably ignore since from the
        data structure point of view it is not really a problem that a vertex is
        incident on two holes - a hole can be seen as a special type of face.
        The function returns false if the index of the face is not valid,
        otherwise the function must complete. """
        return lib_py_gel.Manifold_remove_face(self.obj, fid)
    def remove_edge(self,hid):
        """ Remove an edge from the Manifold. This function will remove the
        faces on either side and the edge itself in the process. Thus, it is a
        simple application of remove_face. """
        return lib_py_gel.Manifold_remove_edge(self.obj, hid)
    def vertex_in_use(self,vid):
        """ check if vertex is in use. This function returns true if the id corresponds
        to a vertex that is currently in the mesh and false otherwise. The id could
        be outside the range of used ids and it could also correspond to a vertex
        which is not active. The function returns false in both cases. """
        return lib_py_gel.Manifold_vertex_in_use(self.obj, vid)
    def face_in_use(self,fid):
        """ check if face is in use. This function returns true if the id corresponds
        to a face that is currently in the mesh and false otherwise. The id could
        be outside the range of used ids and it could also correspond to a face
        which is not active. The function returns false in both cases. """
        return lib_py_gel.Manifold_face_in_use(self.obj, fid)
    def halfedge_in_use(self,hid):
        """ check if halfedge is in use. This function returns true if the id corresponds
        to a halfedge that is currently in the mesh and false otherwise. The id could
        be outside the range of used ids and it could also correspond to a halfedge
        which is not active. The function returns false in both cases. """
        return lib_py_gel.Manifold_halfedge_in_use(self.obj, hid)
    def flip_edge(self,hid):
        """ Flip the edge separating two faces. The function first verifies that
        the edge is flippable. This entails making sure that all of the
        following are true.
        1.  adjacent faces are triangles.
        2. neither end point has valency three or less.
        3. the vertices that will be connected are not already.
        If the tests are passed, the flip is performed and the function
        returns True. Otherwise False."""
        return lib_py_gel.Manifold_flip_edge(self.obj,hid)
    def collapse_edge(self,hid, avg_vertices=False):
        """ Collapse an edge.
        Before collapsing, a number of tests are made:
        ---
        1.  For the two vertices adjacent to the edge, we generate a list of all their neighbouring vertices.
        We then generate a  list of the vertices that occur in both these lists.
        That is, we find all vertices connected by edges to both endpoints of the edge and store these in a list.
        2.  For both faces incident on the edge, check whether they are triangular.
        If this is the case, the face will be removed, and it is ok that the the third vertex is connected to both endpoints.
        Thus the third vertex in such a face is removed from the list generated in 1.
        3.  If the list is now empty, all is well.
        Otherwise, there would be a vertex in the new mesh with two edges connecting it to the same vertex. Return false.
        4.  TETRAHEDRON TEST:
        If the valency of both vertices is three, and the incident faces are triangles, we also disallow the operation.
        Reason: A vertex valency of two and two triangles incident on the adjacent vertices makes the construction collapse.
        5.  VALENCY 4 TEST:
        If a triangle is adjacent to the edge being collapsed, it disappears.
        This means the valency of the remaining edge vertex is decreased by one.
        A valency two vertex reduced to a valency one vertex is considered illegal.
        6.  PREVENT MERGING HOLES:
        Collapsing an edge with boundary endpoints and valid faces results in the creation where two holes meet.
        A non manifold situation. We could relax this...
        7. New test: if the same face is in the one-ring of both vertices but not adjacent to the common edge,
        then the result of a collapse would be a one ring where the same face occurs twice. This is disallowed as the resulting
        face would be non-simple.
        If the tests are passed, the collapse is performed and the function
        returns True. Otherwise False."""
        return lib_py_gel.Manifold_collapse_edge(self.obj, hid, avg_vertices)
    def split_face_by_edge(self,fid,v0,v1):
        """   Split a face. The face, f, is split by creating an edge with
        endpoints v0 and v1 (the next two arguments). The vertices of the old
        face between v0 and v1 (in counter clockwise order) continue to belong
        to f. The vertices between v1 and v0 belong to the new face. A handle to
        the new face is returned. """
        return lib_py_gel.Manifold_split_face_by_edge(self.obj, fid, v0, v1)
    def split_face_by_vertex(self,fid):
        """   Split a polygon, f, by inserting a vertex at the barycenter. This
        function is less likely to create flipped triangles than the
        split_face_triangulate function. On the other hand, it introduces more
        vertices and probably makes the triangles more acute. A handle to the
        inserted vertex is returned. """
        return lib_py_gel.Manifold_split_face_by_vertex(self.obj,fid)
    def split_edge(self,hid):
        """   Insert a new vertex on halfedge h. The new halfedge is insterted
        as the previous edge to h. A handle to the inserted vertex is returned. """
        return lib_py_gel.Manifold_split_edge(self.obj,hid)
    def stitch_boundary_edges(self,h0,h1):
        """   Stitch two halfedges. Two boundary halfedges can be stitched
        together. This can be used to build a complex mesh from a bunch of
        simple faces. """
        return lib_py_gel.Manifold_stitch_boundary_edges(self.obj, h0, h1)
    def merge_faces(self,hid):
        """   Merges two faces into a single polygon. The first face is f. The
        second face is adjacent to f along the halfedge h. This function returns
        true if the merging was possible and false otherwise. Currently merge
        only fails if the mesh is already illegal. Thus it should, in fact,
        never fail. """
        if self.is_halfedge_at_boundary(hid):
            return False
        fid = self.incident_face(hid)
        return lib_py_gel.Manifold_merge_faces(self.obj, fid, hid)
    def close_hole(self,hid):
        """ Close hole given by hid (i.e. the face referenced by hid). Returns
        index of the created face or the face that was already there if, in
        fact, hid was not next to a hole. """
        return lib_py_gel.Manifold_close_hole(self.obj, hid)
    def cleanup(self):
        """ Remove unused items from Mesh, map argument is to be used for
        attribute vector cleanups in order to maintain sync."""
        lib_py_gel.Manifold_cleanup(self.obj)
    def is_halfedge_at_boundary(self, hid):
        """ Returns True if halfedge is a boundary halfedge, i.e. face on either
        side is invalid. """
        return lib_py_gel.is_halfedge_at_boundary(self.obj, hid)
    def is_vertex_at_boundary(self, vid):
        """ Returns True if vertex lies at a boundary. """
        return lib_py_gel.is_vertex_at_boundary(self.obj, vid)
    def edge_length(self, hid):
        """ Returns length of edge passed as argument. """
        return lib_py_gel.length(self.obj, hid)
    # def boundary_edge(self,vid):
    #     hid = 0
    #     if not lib_py_gel.boundary_edge(self.obj,vid,hid):
    #         return None
    #     return hid
    def valency(self,vid):
        """ Returns valency, i.e. number of incident edges."""
        return lib_py_gel.valency(self.obj,vid)
    def vertex_normal(self, vid):
        """ Returns vertex normal (angle weighted) """
        n = (ct.c_double*3)()
        lib_py_gel.vertex_normal(self.obj, vid, ct.byref(n))
        return np.array([n[0],n[1],n[2]])
    def connected(self, v0, v1):
        """ Returns true if the two argument vertices are in each other's one-rings."""
        return lib_py_gel.connected(self.obj,v0,v1)
    def no_edges(self, fid):
        """ Compute the number of edges of a face """
        return lib_py_gel.no_edges(self.obj, fid)
    def face_normal(self, fid):
        """ Compute the normal of a face. If the face is not a triangle,
        the normal is not defined, but computed using the first three
        vertices of the face. """
        n = (ct.c_double*3)()
        lib_py_gel.face_normal(self.obj, fid, ct.byref(n))
        return np.array([n[0],n[1],n[2]])
    def area(self, fid):
        """ Returns the area of a face. """
        return lib_py_gel.area(self.obj, fid)
    def perimeter(self, fid):
        """ Returns the perimeter of a face. """
        return lib_py_gel.perimeter(self.obj, fid)
    def centre(self, fid):
        """ Returns the centre of a face. """
        v = (ct.c_double*3)()
        lib_py_gel.centre(self.obj, fid, ct.byref(v))
        return v

lib_py_gel.valid.restype = ct.c_bool
lib_py_gel.valid.argtypes = (ct.c_void_p,)
def valid(m):
    """  Verify Manifold Integrity Performs a series of tests to check that this
    is a valid manifold. This function is not rigorously constructed but seems
    to catch all problems so far. The function returns true if the mesh is valid
    and false otherwise. """
    return lib_py_gel.valid(m.obj)

lib_py_gel.closed.restype = ct.c_bool
lib_py_gel.closed.argtypes = (ct.c_void_p,)
def closed(m):
    """ Returns true if the mesh is closed, i.e. has no boundary."""
    return lib_py_gel.closed(m.obj)

lib_py_gel.bbox.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double*3),ct.POINTER(ct.c_double*3))
def bbox(m):
    """ Returns the min and max corners of the bounding box of the manifold. """
    pmin = (ct.c_double*3)()
    pmax = (ct.c_double*3)()
    lib_py_gel.bbox(m.obj, ct.byref(pmin),ct.byref(pmax))
    return (np.ctypeslib.as_array(pmin,3),np.ctypeslib.as_array(pmax,3))

lib_py_gel.bsphere.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double*3), ct.POINTER(ct.c_double))
def bsphere(m):
    """ Calculate the bounding sphere of the manifold.
    Returns centre,radius """
    c = (ct.c_double*3)()
    r = (ct.c_double)()
    lib_py_gel.bsphere(m.obj,ct.byref(c),ct.byref(r))
    return (c,r)

lib_py_gel.stitch_mesh.argtypes = (ct.c_void_p,ct.c_double)
lib_py_gel.stitch_mesh.restype = ct.c_int
def stitch(m, rad=1e-30):
    """ Stitch together edges whose endpoints coincide geometrically. This
    function allows you to create a mesh as a bunch of faces and then stitch
    these together to form a coherent whole. What this function adds is a
    spatial data structure to find out which vertices coincide. The return value
    is the number of edges that could not be stitched. Often this is because it
    would introduce a non-manifold situation."""
    return lib_py_gel.stitch_mesh(m.obj,rad)

lib_py_gel.obj_save.argtypes = (ct.c_char_p, ct.c_void_p)
def obj_save(fn, m):
    """ Save Manifold to Wavefront obj file. """
    s = ct.c_char_p(fn.encode('utf-8'))
    lib_py_gel.obj_save(s, m.obj)

lib_py_gel.off_save.argtypes = (ct.c_char_p, ct.c_void_p)
def off_save(fn, m):
    """ Save Manifold to OFF file. """
    s = ct.c_char_p(fn.encode('utf-8'))
    lib_py_gel.off_save(s, m.obj)

lib_py_gel.x3d_save.argtypes = (ct.c_char_p, ct.c_void_p)
def x3d_save(fn, m):
    """ Save Manifold to X3D file. """
    s = ct.c_char_p(fn.encode('utf-8'))
    lib_py_gel.x3d_save(s, m.obj)

lib_py_gel.obj_load.argtypes = (ct.c_char_p, ct.c_void_p)
def obj_load(fn):
    """ Load Manifold from Wavefront obj file. """
    m = Manifold()
    s = ct.c_char_p(fn.encode('utf-8'))
    if lib_py_gel.obj_load(s, m.obj):
        return m
    return None

lib_py_gel.off_load.argtypes = (ct.c_char_p, ct.c_void_p)
def off_load(fn):
    """ Load Manifold from OFF file. """
    m = Manifold()
    s = ct.c_char_p(fn.encode('utf-8'))
    if lib_py_gel.off_load(s, m.obj):
        return m
    return None

lib_py_gel.ply_load.argtypes = (ct.c_char_p, ct.c_void_p)
def ply_load(fn):
    """ Load Manifold from Stanford PLY file. """
    m = Manifold()
    s = ct.c_char_p(fn.encode('utf-8'))
    if lib_py_gel.ply_load(s, m.obj):
        return m
    return None

lib_py_gel.x3d_load.argtypes = (ct.c_char_p, ct.c_void_p)
def x3d_load(fn):
    """ Load Manifold from X3D file. """
    m = Manifold()
    s = ct.c_char_p(fn.encode('utf-8'))
    if lib_py_gel.x3d_load(s, m.obj):
        return m
    return None

lib_py_gel.remove_caps.argtypes = (ct.c_void_p, ct.c_float)
def remove_caps(m, thresh=2.9):
    """ Remove caps from a manifold consisting of only triangles. A cap is a
    triangle with two very small angles and an angle close to pi, however a cap
    does not necessarily have a very short edge. Set the ang_thresh to a value
    close to pi. The closer to pi the _less_ sensitive the cap removal. A cap is
    removed by flipping the (long) edge E opposite to the vertex V with the
    angle close to pi. However, the function is more complex. Read code and
    document more carefully !!! """
    lib_py_gel.remove_caps(m.obj,thresh)

lib_py_gel.remove_needles.argtypes = (ct.c_void_p, ct.c_float, ct.c_bool)
def remove_needles(m, thresh=0.05, average_positions=False):
    """  Remove needles from a manifold consisting of only triangles. A needle
    is a triangle with a single very short edge. It is moved by collapsing the
    short edge. The thresh parameter sets the length threshold."""
    abs_thresh = thresh * average_edge_length(m)
    lib_py_gel.remove_needles(m.obj,abs_thresh, average_positions)

lib_py_gel.close_holes.argtypes = (ct.c_void_p,ct.c_int)
def close_holes(m, max_size=100):
    """  This function replaces holes by faces. It is really a simple function
    that just finds all loops of edges next to missing faces. """
    lib_py_gel.close_holes(m.obj, max_size)

lib_py_gel.flip_orientation.argtypes = (ct.c_void_p,)
def flip_orientation(m):
    """  Flip the orientation of a mesh. After calling this function, normals
    will point the other way and clockwise becomes counter clockwise """
    lib_py_gel.flip_orientation(m.obj)

lib_py_gel.merge_coincident_boundary_vertices.argtypes = (ct.c_void_p, ct.c_double)
def merge_coincident_boundary_vertices(m, rad = 1.0e-30):
    """  Merg vertices that are boundary vertices and coincident.
        However, if one belongs to the other's one ring or the onr
        rings share a vertex, they will not be merged. """
    lib_py_gel.merge_coincident_boundary_vertices(m.obj, rad)

lib_py_gel.minimize_curvature.argtypes = (ct.c_void_p,ct.c_bool)
def minimize_curvature(m,anneal=False):
    """ Minimizes mean curvature. This is really the same as dihedral angle
    minimization, except that we weight by edge length. """
    lib_py_gel.minimize_curvature(m.obj, anneal)

lib_py_gel.minimize_dihedral_angle.argtypes = (ct.c_void_p, ct.c_int, ct.c_bool, ct.c_bool, ct.c_double)
def minimize_dihedral_angle(m,max_iter=10000, anneal=False, alpha=False, gamma=4.0):
    """ Minimizes dihedral angles.
        Arguments:
        max_iter is the maximum number of iterations for simulated annealing.
        anneal tells us the code whether to apply simulated annealing
        alpha=False means that we use the cosine of angles rather than true angles (faster)
        gamma is the power to which the angles are raised."""
    lib_py_gel.minimize_dihedral_angle(m.obj, max_iter, anneal,alpha,ct.c_double(gamma))


lib_py_gel.maximize_min_angle.argtypes = (ct.c_void_p,ct.c_float,ct.c_bool)
def maximize_min_angle(m,dihedral_thresh=0.95,anneal=False):
    """ Maximizes the minimum angle of triangles. Makes the mesh more Delaunay."""
    lib_py_gel.maximize_min_angle(m.obj,dihedral_thresh,anneal)

lib_py_gel.optimize_valency.argtypes = (ct.c_void_p,ct.c_bool)
def optimize_valency(m,anneal=False):
    """ Tries to achieve valence 6 internally and 4 along edges. """
    lib_py_gel.optimize_valency(m.obj, anneal)

lib_py_gel.randomize_mesh.argtypes = (ct.c_void_p,ct.c_int)
def randomize_mesh(m,max_iter=1):
    """  Make random flips. Useful for generating synthetic test cases. """
    lib_py_gel.randomize_mesh(m.obj, max_iter)

lib_py_gel.quadric_simplify.argtypes = (ct.c_void_p,ct.c_double,ct.c_double,ct.c_bool)
def quadric_simplify(m,keep_fraction,singular_thresh=1e-4,optimal_positions=True):
    """ Garland Heckbert simplification in our own implementation. keep_fraction
    is the fraction of vertices to retain. The singular_thresh defines how small
    singular values from the SVD we accept. It is relative to the greatest
    singular value. If optimal_positions is true, we reposition vertices.
    Otherwise the vertices are a subset of the old vertices."""
    lib_py_gel.quadric_simplify(m.obj, keep_fraction, singular_thresh,optimal_positions)

lib_py_gel.average_edge_length.argtypes = (ct.c_void_p,)
lib_py_gel.average_edge_length.restype = ct.c_float
def average_edge_length(m,max_iter=1):
    """ Returns the average edge length. """
    return lib_py_gel.average_edge_length(m.obj)

lib_py_gel.median_edge_length.argtypes = (ct.c_void_p,)
lib_py_gel.median_edge_length.restype = ct.c_float
def median_edge_length(m,max_iter=1):
    """ Returns the median edge length """
    return lib_py_gel.median_edge_length(m.obj)

lib_py_gel.refine_edges.argtypes = (ct.c_void_p,ct.c_float)
lib_py_gel.refine_edges.restype = ct.c_int
def refine_edges(m,threshold):
    """ Split all edges in mesh passed as first argument which are longer
    than the threshold (second arg) length. A split edge
    results in a new vertex of valence two."""
    return lib_py_gel.refine_edges(m.obj, threshold)

lib_py_gel.cc_split.argtypes = (ct.c_void_p,)
def cc_split(m):
    """ Perform a Catmull-Clark split, i.e. a split where each face is divided
    into new quadrilateral faces formed by connecting a corner with a point on
    each incident edge and a point at the centre of the face."""
    lib_py_gel.cc_split(m.obj)

lib_py_gel.loop_split.argtypes = (ct.c_void_p,)
def loop_split(m):
    """ Perform a loop split where each edge is divided into two segments, and
    four new triangles are created for each original triangle. """
    lib_py_gel.loop_split(m.obj)

lib_py_gel.root3_subdivide.argtypes = (ct.c_void_p,)
def root3_subdivide(m):
    """ Leif Kobbelt's subdivision scheme, where a vertex is positions in the
    center of each face and all old edges are flipped. """
    lib_py_gel.root3_subdivide(m.obj)

lib_py_gel.rootCC_subdivide.argtypes = (ct.c_void_p,)
def rootCC_subdivide(m):
    """ This subd. scheme creates a vertex inside each original (quad) face,
    producing four triangles. Triangles sharing an old edge are then merged.
    Two steps produce something similar to Catmull-Clark. """
    lib_py_gel.rootCC_subdivide(m.obj)

lib_py_gel.butterfly_subdivide.argtypes = (ct.c_void_p,)
def butterfly_subdivide(m):
    """ An interpolatory scheme. Creates the same connectivity as Loop. """
    lib_py_gel.butterfly_subdivide(m.obj)

lib_py_gel.cc_smooth.argtypes = (ct.c_void_p,)
def cc_smooth(m):
    """ If called after cc_split, this function completes a Catmull-Clark
    subdivision step. """
    lib_py_gel.cc_smooth(m.obj)

lib_py_gel.loop_smooth.argtypes = (ct.c_void_p,)
def loop_smooth(m):
    """ If called after Loop split, this function completes a step of Loop
    subdivision. """
    lib_py_gel.loop_smooth(m.obj)

lib_py_gel.ear_clip_triangulate.argtypes = (ct.c_void_p,)
lib_py_gel.shortest_edge_triangulate.argtypes = (ct.c_void_p,)
def triangulate(m, clip_ear=True):
    """ Turn a general polygonal mesh into a triangle mesh by repeatedly
        splitting a polygon into smaller polygons. """
    if clip_ear:
        lib_py_gel.ear_clip_triangulate(m.obj)
    else:
        lib_py_gel.shortest_edge_triangulate(m.obj)


try:
    lib_py_gel.GLManifoldViewer_new.restype = ct.c_void_p
    lib_py_gel.GLManifoldViewer_delete.argtypes = (ct.c_void_p,)
    lib_py_gel.GLManifoldViewer_display.argtypes = (ct.c_void_p,ct.c_void_p,ct.c_void_p,ct.c_char,ct.c_bool, ct.POINTER(ct.c_float*3), ct.POINTER(ct.c_double),ct.c_bool,ct.c_bool)
    lib_py_gel.GLManifoldViewer_get_annotation_points.restype = ct.c_size_t
    lib_py_gel.GLManifoldViewer_get_annotation_points.argtypes = (ct.c_void_p, ct.POINTER(ct.POINTER(ct.c_double)))
    lib_py_gel.GLManifoldViewer_set_annotation_points.argtypes = (ct.c_void_p, ct.c_int, ct.POINTER(ct.c_double))
    lib_py_gel.GLManifoldViewer_event_loop.argtypes = (ct.c_bool,)
    class Viewer:
        """ An OpenGL Viewer for Manifolds and Graphs. Having created an instance of this
        class, call display to show a mesh or a graph. The display function is flexible,
        allowing several types of interactive visualization. Each instance of this
        class corresponds to a single window, but you can have several
        GLManifoldViewer and hence also several windows showing different
        visualizations. """
        def __init__(self):
            self.obj = lib_py_gel.GLManifoldViewer_new()
        def __del__(self):
            lib_py_gel.GLManifoldViewer_delete(self.obj)
        def display(self, m, g=None, mode='w', smooth=True, bg_col=[0.3,0.3,0.3], data=None, reset_view=False, once=False):
            """ Display a mesh

            Args:
            ---
            - m : the Manifold mesh or Graph we want to show.
            - g : the Graph we want to show. If you only want to show a graph, you
                can simply pass the graph as m, so the g argument is relevant only if
                you need to show both a Manifold _and_ a Graph.
            - mode : a single character that determines how the mesh is visualized:
                'w' - wireframe,
                'i' - isophote,
                'g' - glazed (try it and see),
                's' - scalar field,
                'l' - line field,
                'n' - normal.
                'x' - xray or ghost rendering. Useful to show Manifold on top of Graph
            - smooth : if True we use vertex normals. Otherwise, face normals.
            - bg_col : background color.
            - data : per vertex data for visualization. scalar or vector field.
            - reset_view : if False view is as left in the previous display call. If
                True, the view is reset to the default.
            - once : if True we immediately exit the event loop and return. However,
                the window stays and if the event loop is called from this or any
                other viewer, the window will still be responsive.

            Interactive controls:
            ---
            When a viewer window is displayed on the screen, you can naviagate with
            the mouse: Left mouse button rotates, right mouse button is used for
            zooming and (if shift is pressed) for panning. If you hold control, any
            mouse button will pick a point on the 3D model. Up to 19 of these points
            have unique colors.  If you pick an already placed annotation point it
            will be removed and can now be placed elsewhere. Hit space bar to clear
            the annotation points. Hitting ESC exits the event loop causing control
            to return to the script.
            """
            data_ct = np.array(data,dtype=ct.c_double).ctypes
            data_a = data_ct.data_as(ct.POINTER(ct.c_double))
            bg_col_ct = np.array(bg_col,dtype=ct.c_float).ctypes
            bg_col_a = bg_col_ct.data_as(ct.POINTER(ct.c_float*3))
            if isinstance(m,Graph):
                g = m
                m = None
            if isinstance(m,Manifold) and isinstance(g,Graph):
                lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, g.obj, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
            elif isinstance(m,Manifold):
                lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, 0, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
            elif isinstance(g,Graph):
                lib_py_gel.GLManifoldViewer_display(self.obj, 0, g.obj, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
                
        def annotation_points(self):
            """ Retrieve a vector of annotation points. This vector is not a copy,
            so any changes made to the points will be reflected in the viewer. """
            pos = ct.POINTER(ct.c_double)()
            n = lib_py_gel.GLManifoldViewer_get_annotation_points(self.obj, ct.byref(pos))
            if n == 0:
                return None
            return np.ctypeslib.as_array(pos,(n,3))
        def set_annotation_points(self, pts):
            n = int(np.size(pts)/3)
            pts_ct = np.array(pts,dtype=ct.c_double).ctypes
            pts_a = pts_ct.data_as(ct.POINTER(ct.c_double))
            lib_py_gel.GLManifoldViewer_set_annotation_points(self.obj, n, pts_a)
        @staticmethod
        def event_loop():
            """ Explicit call to the event loop. This function enters the event loop.
            Call it if you want to turn on interactivity in the currently displayed
            window."""
            lib_py_gel.GLManifoldViewer_event_loop(False)
except AttributeError:
    pass

lib_py_gel.MeshDistance_new.restype = ct.c_void_p
lib_py_gel.MeshDistance_new.argtypes = (ct.c_void_p,)
lib_py_gel.MeshDistance_signed_distance.argtypes = (ct.c_void_p,ct.c_int, ct.POINTER(ct.c_float),ct.POINTER(ct.c_float),ct.c_float)
lib_py_gel.MeshDistance_ray_inside_test.argtypes = (ct.c_void_p,ct.c_int, ct.POINTER(ct.c_float),ct.POINTER(ct.c_int),ct.c_int)

lib_py_gel.MeshDistance_delete.argtypes = (ct.c_void_p,)
class MeshDistance:
    """ This class allows you to compute the distance from any point in space to
    a Manifold. The constructor creates an instance based on a specific mesh,
    and the signed_distance function computes the actual distance. """
    def __init__(self,m):
        self.obj = lib_py_gel.MeshDistance_new(m.obj)
    def __del__(self):
        lib_py_gel.MeshDistance_delete(self.obj)
    def signed_distance(self,pts,upper=1e30):
        """ Compute the signed distance from p to the mesh stored in this class
        instance. The distance is positive if outside and negative inside. The
        upper parameter can be used to threshold how far away the distance is of
        interest. """
        p = np.reshape(np.array(pts,dtype=ct.c_float), (-1,3))
        n = p.shape[0]
        d = np.ndarray(n, dtype=ct.c_float)
        p_ct = p.ctypes.data_as(ct.POINTER(ct.c_float))
        d_ct = d.ctypes.data_as(ct.POINTER(ct.c_float))
        lib_py_gel.MeshDistance_signed_distance(self.obj,n,p_ct,d_ct,upper)
        return d
    def ray_inside_test(self,pts,no_rays=3):
        """Check whether a point is inside or outside the stored by casting rays.
        Effectively, this is the sign of the distance. In some cases casting (multiple)
        ray is more robust than using the sign computed locally. """
        p = np.reshape(np.array(pts,dtype=ct.c_float), (-1,3))
        n = p.shape[0]
        s = np.ndarray(n, dtype=ct.c_int)
        p_ct = p.ctypes.data_as(ct.POINTER(ct.c_float))
        s_ct = s.ctypes.data_as(ct.POINTER(ct.c_int))
        lib_py_gel.MeshDistance_ray_inside_test(self.obj,n,p_ct,s_ct,no_rays)
        return s


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

class Graph:
    """ This class is for representing graphs embedded in 3D. The class does not in
    itself come with many features: it contains methods for creating, accessing, and
    housekeeping. """
    def __init__(self,orig=None):
        if orig == None:
            self.obj = lib_py_gel.Graph_new()
        else:
            self.obj = lib_py_gel.Graph_copy(orig.obj)
    def __del__(self):
        lib_py_gel.Graph_delete(self.obj)
    def clear(self):
        """ Clear the graph. """
        lib_py_gel.Graph_clear(self.obj)
    def cleanup(self):
        """ Cleanup reorders the graph nodes such that there is no
        gap in the index range. """
        lib_py_gel.Graph_cleanup(self.obj)
    def nodes(self):
        """ Get all nodes as an iterable range """
        nodes = IntVector()
        lib_py_gel.Graph_nodes(self.obj, nodes.obj)
        return nodes
    def neighbors(self, n, mode='n'):
        """ Get the neighbors of node n. The final argument is either 'n' or 'e'. If it is 'n'
        the function returns all neighboring nodes, and if it is 'e' it returns incident edges."""
        nbors = IntVector()
        lib_py_gel.Graph_neighbors(self.obj, n, nbors.obj, ct.c_char(mode.encode('ascii')))
        return nbors
    def positions(self):
        """ Get the vertex positions by reference. You can assign to the
        positions. """
        pos = ct.POINTER(ct.c_double)()
        n = lib_py_gel.Graph_positions(self.obj, ct.byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
    def average_edge_length(self):
        """ Returns the average edge length. """
        ael = lib_py_gel.Graph_average_edge_length(self.obj)
        return ael
    def add_node(self, p):
        """ Adds node with position p to the graph and returns the
        index of the new node. """
        return lib_py_gel.Graph_add_node(self.obj, np.array(p))
    def remove_node(self, n):
        """ Removes the node n passed as argument. This does not change
        any indices of other nodes, but n is then invalid. """
        lib_py_gel.Graph_remove_node(self.obj, n)
    def node_in_use(self, n):
        """ Checks if n is in_use. This function returns false both
        if n has been removed and if n is an index outside the range of
        indices that are used. """
        return lib_py_gel.Graph_node_in_use(self.obj, n)
    def connect_nodes(self, n0, n1):
        """ Creates a new edge connecting nodes n0 and n1. The index of
        the new edge is returned. """
        return lib_py_gel.Graph_connect_nodes(self.obj, n0, n1)
    def disconnect_nodes(self, n0, n1):
        """ Disconect nodes n0 and n1"""
        lib_py_gel.Graph_disconnect_nodes(self.obj, n0, n1)
    def merge_nodes(self, n0, n1, avg_pos):
        """ Merge nodes n0 and n1. avg_pos indicates if you want the position to be the average. """
        lib_py_gel.Graph_merge_nodes(self.obj, n0, n1, avg_pos)


lib_py_gel.graph_from_mesh.argtypes = (ct.c_void_p, ct.c_void_p)
def graph_from_mesh(m):
    """ Creates a graph from a mesh. The argument, m, is the input mesh,
    and the function returns a graph with the same vertices and edges
    as m."""
    g = Graph()
    lib_py_gel.graph_from_mesh(m.obj, g.obj)
    return g

lib_py_gel.graph_load.argtypes = (ct.c_void_p, ct.c_char_p)
lib_py_gel.graph_load.restype = ct.c_void_p
def graph_load(fn):
    """ Load a graph from a file. The argument, fn, is the filename which
    is in a special format similar to Wavefront obj. The loaded graph is
    returned by the function - or None if loading failed. """
    s = ct.c_char_p(fn.encode('utf-8'))
    g = Graph();
    if lib_py_gel.graph_load(g.obj, s):
        return g
    return None

lib_py_gel.graph_save.argtypes = (ct.c_void_p, ct.c_char_p)
lib_py_gel.graph_save.restype = ct.c_bool
def graph_save(fn, g):
    """ Save graph to a file. The first argument, fn, is the file name,
    and g is the graph. This function returns True if saving happened and
    False otherwise. """
    s = ct.c_char_p(fn.encode('utf-8'))
    return lib_py_gel.graph_save(g.obj, s)

lib_py_gel.graph_to_mesh_cyl.argtypes = (ct.c_void_p, ct.c_void_p, ct.c_float)
lib_py_gel.graph_to_mesh_cyl.restype = ct.c_void_p
def graph_to_mesh_cyl(g, fudge):
    """ Creates a Manifold mesh from the graph. The first argument, g, is the
    mesh we want converted, and fudge is a constant that is used to increase the radius
    of every node. This is useful if the radii are 0. """
    m = Manifold()
    lib_py_gel.graph_to_mesh_cyl(g.obj, m.obj, fudge)
    return m

lib_py_gel.graph_smooth.argtypes = (ct.c_void_p, ct.c_int, ct.c_float)
def graph_smooth(g, iter=1, alpha=1.0):
    """ Simple Laplacian smoothing of a graph. The first argument is the Graph, g, iter
    is the number of iterations, and alpha is the weight. If the weight is high,
    each iteration causes a lot of smoothing, and a high number of iterations
    ensures that the effect of smoothing diffuses throughout the graph, i.e. that the
    effect is more global than local. """
    lib_py_gel.graph_smooth(g.obj, iter, alpha)

lib_py_gel.graph_edge_contract.argtypes = (ct.c_void_p, ct.c_double)
def graph_edge_contract(g, dist_thresh):
    """ Simplified a graph by contracting edges. The first argument, g, is the graph,
    and only edges shorter than dist_thresh are contracted. When an edge is contracted
    the merged vertices are moved to the average of their former positions. Thus,
    the ordering in which contractions are carried out matters. Hence, edges are contracted
    in the order of increasing length and edges are only considered if neither end point
    is the result of a contraction, but the process is then repeated until no more contractions
    are possible. Returns total number of contractions. """
    return lib_py_gel.graph_edge_contract(g.obj, dist_thresh)

lib_py_gel.graph_prune.argtypes = (ct.c_void_p,)
def graph_prune(g):
    """ Prune leaves of a graph. The graph, g, is passed as the argument. This function removes
    leaf nodes (valency 1) whose only neighbour has valency > 2. In practice such isolated
    leaves are frequently spurious if the graph is a skeleton. Does not return a value. """
    lib_py_gel.graph_prune(g.obj)
    
lib_py_gel.graph_LS_skeleton.argtypes = (ct.c_void_p, ct.c_void_p, ct.c_bool)
def graph_LS_skeleton(g, sampling=True):
    """ Skeletonize a graph using the local separators approach. The first argument, g, is
    the graph, and, sampling indicates whether we try to use all vertices (False) as starting
    points for finding separators or just a sampling (True). The function returns a new graph
    which is the skeleton of the input graph. """
    skel = Graph()
    lib_py_gel.graph_LS_skeleton(g.obj, skel.obj, sampling)
    return skel
