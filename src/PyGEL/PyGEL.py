from ctypes import *
#from numpy.ctypeslib import as_array
import numpy as np
import random
from copy import copy
from numpy.linalg import norm
import sys, os
pathname = os.path.dirname(sys.argv[0])
lib = cdll.LoadLibrary("libPyGEL.dylib")

lib.IntVector_new.restype = c_void_p
lib.IntVector_get.argtypes = (c_void_p, c_size_t)
lib.IntVector_size.argtypes = (c_void_p,)
lib.IntVector_size.restype = c_size_t
lib.IntVector_delete.argtypes = (c_void_p,)
class IntVector:
    def __init__(self):
        self.obj = lib.IntVector_new(0)
    def __del__(self):
        lib.IntVector_delete(self.obj)
    def __len__(self):
        return int(lib.IntVector_size(self.obj))
    def __getitem__(self,key):
        return lib.IntVector_get(self.obj,key)
    def __iter__(self):
        n = lib.IntVector_size(self.obj)
        for i in range(0,n):
            yield lib.IntVector_get(self.obj, i)

lib.Vec3dVector_new.restype = c_void_p
lib.Vec3dVector_get.argtypes = (c_void_p, c_size_t)
lib.Vec3dVector_get.restype = POINTER(c_double)
lib.Vec3dVector_size.argtypes = (c_void_p,)
lib.Vec3dVector_size.restype = c_size_t
lib.Vec3dVector_delete.argtypes = (c_void_p,)
class Vec3dVector:
    def __init__(self):
        self.obj = lib.Vec3dVector_new(0)
    def __del__(self):
        lib.Vec3dVector_delete(self.obj)
    def __len__(self):
        return int(lib.Vec3dVector_size(self.obj))
    def __getitem__(self,key):
        return lib.Vec3dVector_get(self.obj,key)
    def __iter__(self):
        n = lib.Vec3dVector_size(self.obj)
        for i in range(0,n):
            data = lib.Vec3dVector_get(self.obj, i)
            yield [data[0], data[1], data[2]]


lib.I3DTree_new.restype = c_void_p
lib.I3DTree_delete.argtypes = (c_void_p,)
lib.I3DTree_insert.argtypes = (c_void_p, c_double, c_double, c_double, c_size_t)
lib.I3DTree_build.argtypes = (c_void_p,)
lib.I3DTree_closest_point.argtypes = (c_void_p, c_double, c_double, c_double, c_double, POINTER(c_double*3), POINTER(c_size_t))
lib.I3DTree_in_sphere.argtypes = (c_void_p, c_double, c_double, c_double, c_double, c_void_p,c_void_p)
class I3DTree:
    def __init__(self):
        self.obj = lib.I3DTree_new()
    def __del__(self):
        lib.I3DTree_delete(self.obj)
    def insert(self,p,v):
        lib.I3DTree_insert(self.obj, p[0],p[1],p[2],v)
    def build(self):
        lib.I3DTree_build(self.obj)
    def closest_point(self,p,r):
        key = (c_double * 3)()
        val = c_size_t()
        n = lib.I3DTree_closest_point(self.obj, p[0],p[1],p[2],r,byref(key),byref(val))
        if n==1:
            return ([key[0],key[1],key[2]],val.value)
        return None
    def in_sphere(self, p, r):
        keys = Vec3dVector()
        vals = IntVector()
        n = lib.I3DTree_in_sphere(self.obj, p[0],p[1],p[2],r,keys.obj,vals.obj)
        return (keys,vals)


lib.Manifold_new.restype = c_void_p
lib.Manifold_copy.restype = c_void_p
lib.Manifold_copy.argtypes = (c_void_p,)
lib.Manifold_delete.argtypes = (c_void_p,)
lib.Manifold_positions.restype = c_size_t
lib.Manifold_positions.argtypes = (c_void_p, POINTER(POINTER(c_double)))
lib.Manifold_no_allocated_vertices.restype = c_size_t
lib.Manifold_no_allocated_vertices.argtypes = (c_void_p,)
lib.Manifold_no_allocated_faces.restype = c_size_t
lib.Manifold_no_allocated_faces.argtypes = (c_void_p,)
lib.Manifold_no_allocated_halfedges.restype = c_size_t
lib.Manifold_no_allocated_halfedges.argtypes = (c_void_p,)
lib.Manifold_vertices.restype = c_size_t
lib.Manifold_vertices.argtypes = (c_void_p, c_void_p)
lib.Manifold_faces.restype = c_size_t
lib.Manifold_faces.argtypes = (c_void_p, c_void_p)
lib.Manifold_halfedges.restype = c_size_t
lib.Manifold_halfedges.argtypes = (c_void_p,c_void_p)
lib.Manifold_circulate_vertex.restype = c_size_t
lib.Manifold_circulate_vertex.argtypes = (c_void_p, c_size_t, c_char, c_void_p)
lib.Manifold_circulate_face.restype = c_size_t
lib.Manifold_circulate_face.argtypes = (c_void_p, c_size_t, c_char, c_void_p)
lib.Manifold_add_face.argtypes = (c_void_p, c_size_t, np.ctypeslib.ndpointer(c_double))
lib.Manifold_remove_face.restype = c_bool
lib.Manifold_remove_face.argtypes = (c_void_p, c_size_t)
lib.Manifold_remove_edge.restype = c_bool
lib.Manifold_remove_edge.argtypes = (c_void_p, c_size_t)
lib.Manifold_remove_vertex.restype = c_bool
lib.Manifold_remove_vertex.argtypes = (c_void_p, c_size_t)
lib.Manifold_vertex_in_use.restype = c_bool
lib.Manifold_vertex_in_use.argtypes = (c_void_p,c_size_t)
lib.Manifold_face_in_use.restype = c_bool
lib.Manifold_face_in_use.argtypes = (c_void_p, c_size_t)
lib.Manifold_halfedge_in_use.restype = c_bool
lib.Manifold_halfedge_in_use.argtypes = (c_void_p, c_size_t)
lib.Manifold_flip_edge.restype =  c_bool
lib.Manifold_flip_edge.argtypes = (c_void_p, c_size_t)
lib.Manifold_collapse_edge.restype = c_bool
lib.Manifold_collapse_edge.argtypes = (c_void_p,c_size_t,c_bool)
lib.Manifold_split_face_by_edge.restype = c_size_t
lib.Manifold_split_face_by_edge.argtypes = (c_void_p, c_size_t,c_size_t,c_size_t)
lib.Manifold_split_face_by_vertex.restype = c_size_t
lib.Manifold_split_face_by_vertex.argtypes = (c_void_p, c_size_t)
lib.Manifold_split_edge.restype = c_size_t
lib.Manifold_split_edge.argtypes = (c_void_p, c_size_t)
lib.Manifold_stitch_boundary_edges.restype = c_bool
lib.Manifold_stitch_boundary_edges.argtypes = (c_void_p, c_size_t,c_size_t)
lib.Manifold_merge_faces.restype = c_bool
lib.Manifold_merge_faces.argtypes = (c_void_p, c_size_t,c_size_t)
lib.Manifold_close_hole.argtypes = (c_void_p,c_size_t)
lib.Manifold_cleanup.argtypes = (c_void_p,)

lib.Walker_next_halfedge.restype = c_size_t
lib.Walker_next_halfedge.argtypes = (c_void_p, c_size_t)
lib.Walker_opposite_halfedge.restype = c_size_t
lib.Walker_opposite_halfedge.argtypes = (c_void_p, c_size_t)
lib.Walker_incident_face.restype = c_size_t
lib.Walker_incident_face.argtypes = (c_void_p, c_size_t)
lib.Walker_incident_vertex.restype = c_size_t
lib.Walker_incident_vertex.argtypes = (c_void_p, c_size_t)

lib.is_halfedge_at_boundary.restype = c_bool
lib.is_halfedge_at_boundary.argtypes = (c_void_p, c_size_t)
lib.is_vertex_at_boundary.restype = c_bool
lib.is_vertex_at_boundary.argtypes = (c_void_p, c_size_t)
lib.length.restype = c_double
lib.length.argtypes = (c_void_p, c_size_t)
lib.boundary_edge.restype = c_bool
lib.boundary_edge.argtypes = (c_void_p, c_size_t, c_size_t)
lib.valency.restype = c_size_t
lib.valency.argtypes = (c_void_p, c_size_t)
lib.vertex_normal.argtypes = (c_void_p, c_size_t, POINTER(c_double*3))
lib.connected.restype = c_bool
lib.connected.argtypes = (c_void_p, c_size_t, c_size_t)
lib.no_edges.restype = c_size_t
lib.no_edges.argtypes = (c_void_p, c_size_t)
lib.face_normal.argtypes = (c_void_p, c_size_t, POINTER(c_double*3))
lib.area.restype = c_double
lib.area.argtypes = (c_void_p, c_size_t)
lib.perimeter.restype = c_double
lib.perimeter.argtypes = (c_void_p, c_size_t)
lib.centre.argtypes = (c_void_p, c_size_t, POINTER(c_double*3))

class Manifold:
    """ The Manifold class represents a halfedge based mesh.
    Since meshes based on the halfedge representation must be manifold (although
    exceptions could be made) the class is thus named. Manifold contains many
    functions for mesh manipulation and associated the position attribute
    with vertices.
    """
    def __init__(self,orig=None):
        if orig == None:
            self.obj = lib.Manifold_new()
        else:
            self.obj = lib.Manifold_copy(orig.obj)
    def __del__(self):
        lib.Manifold_delete(self.obj)
    def add_face(self,pts):
        """ Add a face to the Manifold.
        This function takes a list of points as argument and creates a face
        in the mesh with those points as vertices.
        """
        lib.Manifold_add_face(self.obj, len(pts), np.array(pts))
    def positions(self):
        pos = POINTER(c_double)()
        n = lib.Manifold_positions(self.obj, byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
    def no_allocated_vertices(self):
        return lib.Manifold_no_allocated_vertices(self.obj)
    def no_allocated_faces(self):
        return lib.Manifold_no_allocated_faces(self.obj)
    def no_allocated_halfedges(self):
        return lib.Manifold_no_allocated_halfedges(self.obj)
    def vertices(self):
        verts = IntVector()
        n = lib.Manifold_vertices(self.obj, verts.obj)
        return verts
    def faces(self):
        faces = IntVector()
        n = lib.Manifold_faces(self.obj, faces.obj)
        return faces
    def halfedges(self):
        hedges = IntVector()
        n = lib.Manifold_halfedges(self.obj, hedges.obj)
        return hedges
    def circulate_vertex(self,v,mode='v'):
        nbrs = IntVector()
        n = lib.Manifold_circulate_vertex(self.obj, v, c_char(mode.encode('ascii')), nbrs.obj)
        return nbrs
    def circulate_face(self,f,mode='v'):
        nbrs = IntVector()
        n = lib.Manifold_circulate_face(self.obj, f, c_char(mode.encode('ascii')), nbrs.obj)
        return nbrs
    def next_halfedge(self,hid):
        return lib.Walker_next_halfedge(self.obj, hid)
    def opposite_halfedge(self,hid):
        return lib.Walker_opposite_halfedge(self.obj, hid)
    def incident_face(self,hid):
        return lib.Walker_incident_face(self.obj, hid)
    def incident_vertex(self,hid):
        return lib.Walker_incident_vertex(self.obj, hid)
    def remove_vertex(self,vid):
        return lib.Manifold_remove_vertex(self.obj, vid)
    def remove_face(self,fid):
        return lib.Manifold_remove_face(self.obj, fid)
    def remove_edge(self,hid):
        return lib.Manifold_remove_edge(self.obj, hid)
    def vertex_in_use(self,vid):
        return lib.Manifold_vertex_in_use(self.obj, vid)
    def face_in_use(self,fid):
        return lib.Manifold_face_in_use(self.obj, fid)
    def halfedge_in_use(self,hid):
        return lib.Manifold_halfedge_in_use(self.obj, hid)
    def flip_edge(self,hid):
        return lib.Manifold_flip_edge(self.obj,hid)
    def collapse_edge(self,hid, avg_vertices=False):
        return lib.Manifold_collapse_edge(self.obj, hid, avg_vertices)
    def split_face_by_edge(self,fid,v0,v1):
        return lib.Manifold_split_face_by_edge(self.obj, fid, v0, v1)
    def split_face_by_vertex(self,fid):
        return lib.Manifold_split_face_by_vertex(self.obj,fid)
    def split_edge(self,hid):
        return lib.Manifold_split_edge(self.obj,hid)
    def stitch_boundary_edges(self,h0,h1):
        return lib.Manifold_stitch_boundary_edges(self.obj, h0, h1)
    def merge_faces(self,hid):
        if self.is_halfedge_at_boundary(hid):
            return False
        fid = self.incident_face(hid)
        return lib.Manifold_merge_faces(self.obj, fid, hid)
    def close_hole(self,hid):
        return lib.Manifold_close_hole(self.obj, hid)
    def cleanup(self):
        lib.Manifold_cleanup(self.obj)
    def is_halfedge_at_boundary(self, hid):
        return lib.is_halfedge_at_boundary(self.obj, hid)
    def is_vertex_at_boundary(self, vid):
        return lib.is_vertex_at_boundary(self.obj, vid)
    def edge_length(self, hid):
        return lib.length(self.obj, hid)
    # def boundary_edge(self,vid):
    #     hid = 0
    #     if not lib.boundary_edge(self.obj,vid,hid):
    #         return None
    #     return hid
    def valency(self,vid):
        return lib.valency(self.obj,vid)
    def vertex_normal(self, vid):
        n = (c_double*3)()
        return lib.vertex_normal(self.obj, vid, byref(n))
    def connected(self, v0, v1):
        return lib.connected(self.obj,v0,v1)
    def no_edges(self, fid):
        return lib.no_edges(self.obj, fid)
    def face_normal(self, fid):
        n = (c_double*3)()
        return lib.face_normal(self.obj, fid, byref(n))
    def area(self, fid):
        return lib.area(self.obj, fid)
    def perimeter(self, fid):
        return lib.perimeter(self.obj, fid)
    def centre(self, fid):
        return lib.centre(self, fid)

lib.valid.restype = c_bool
lib.valid.argtypes = (c_void_p,)
def valid(m):
    return lib.valid(m.obj)

lib.closed.restype = c_bool
lib.closed.argtypes = (c_void_p,)
def closed(m):
    return lib.closed(m.obj)

lib.bbox.argtypes = (c_void_p, POINTER(c_double*3),POINTER(c_double*3))
def bbox(m):
    pmin = (c_double*3)()
    pmax = (c_double*3)()
    lib.bbox(m.obj, byref(pmin),byref(pmax))
    return (pmin,pmax)

lib.bsphere.argtypes = (c_void_p, POINTER(c_double*3), POINTER(c_double))
def bsphere(m):
    c = (c_double*3)()
    r = (c_double)()
    lib.bsphere(m.obj,byref(c),byref(r))
    return (c,r)

lib.stitch_mesh.argtypes = (c_void_p,c_double)
def stitch(m, rad=1e-30):
    lib.stitch_mesh(m.obj,rad)

lib.obj_load.argtypes = (c_char_p, c_void_p)
def obj_load(fn):
    m = Manifold()
    s = c_char_p(fn.encode('utf-8'))
    lib.obj_load(s, m.obj)
    return m

lib.off_load.argtypes = (c_char_p, c_void_p)
def off_load(fn):
    m = Manifold()
    s = c_char_p(fn.encode('utf-8'))
    lib.off_load(s, m.obj)
    return m

lib.ply_load.argtypes = (c_char_p, c_void_p)
def ply_load(fn):
    m = Manifold()
    s = c_char_p(fn.encode('utf-8'))
    lib.ply_load(s, m.obj)
    return m

lib.x3d_load.argtypes = (c_char_p, c_void_p)
def x3d_load(fn):
    m = Manifold()
    s = c_char_p(fn.encode('utf-8'))
    lib.x3d_load(s, m.obj)
    return m

lib.remove_caps.argtypes = (c_void_p, c_float)
def remove_caps(m, thresh=2.9):
    lib.remove_caps(m.obj,thresh)

lib.remove_needles.argtypes = (c_void_p, c_float, c_bool)
def remove_needles(m, thresh=0.05, average_positions=False):
    abs_thresh = thresh * average_edge_length(m)
    lib.remove_needles(m.obj,abs_thresh, average_positions)

lib.close_holes.argtypes = (c_void_p,)
def close_holes(m):
    lib.close_holes(m.obj)

lib.flip_orientation.argtypes = (c_void_p,)
def flip_orientation(m):
    lib.flip_orientation(m.obj)

lib.minimize_curvature.argtypes = (c_void_p,c_bool)
def minimize_curvature(m,anneal=False):
    lib.minimize_curvature(m.obj, anneal)

lib.maximize_min_angle.argtypes = (c_void_p,c_float,c_bool)
def maximize_min_angle(m,dihedral_thresh=0.95,anneal=False):
    lib.maximize_min_angle(m.obj,dihedral_thresh,anneal)

lib.optimize_valency.argtypes = (c_void_p,c_bool)
def optimize_valency(m,anneal=False):
    lib.optimize_valency(m.obj, anneal)

lib.randomize_mesh.argtypes = (c_void_p,c_int)
def randomize_mesh(m,max_iter=1):
    lib.randomize_mesh(m.obj, max_iter)

lib.quadric_simplify.argtypes = (c_void_p,c_double,c_double,c_bool)
def quadric_simplify(m,keep_fraction,singular_thresh=1e-4,optimal_positions=True):
    lib.quadric_simplify(m.obj, keep_fraction, singular_thresh,optimal_positions)

lib.average_edge_length.argtypes = (c_void_p,)
lib.average_edge_length.restype = c_float
def average_edge_length(m,max_iter=1):
    return lib.average_edge_length(m.obj)

lib.median_edge_length.argtypes = (c_void_p,)
lib.median_edge_length.restype = c_float
def median_edge_length(m,max_iter=1):
    return lib.median_edge_length(m.obj)

lib.refine_edges.argtypes = (c_void_p,c_float)
lib.refine_edges.restype = c_int
def refine_edges(m,threshold):
    return lib.refine_edges(m.obj, threshold)

lib.cc_split.argtypes = (c_void_p,)
def cc_split(m):
    lib.cc_split(m.obj)

lib.loop_split.argtypes = (c_void_p,)
def loop_split(m):
    lib.loop_split(m.obj)

lib.root3_subdivide.argtypes = (c_void_p,)
def root3_subdivide(m):
    lib.root3_subdivide(m.obj)

lib.rootCC_subdivide.argtypes = (c_void_p,)
def rootCC_subdivide(m):
    lib.rootCC_subdivide(m.obj)

lib.butterfly_subdivide.argtypes = (c_void_p,)
def butterfly_subdivide(m):
    lib.butterfly_subdivide(m.obj)

lib.cc_smooth.argtypes = (c_void_p,)
def cc_smooth(m):
    lib.cc_smooth(m.obj)

lib.loop_smooth.argtypes = (c_void_p,)
def loop_smooth(m):
    lib.loop_smooth(m.obj)

lib.shortest_edge_triangulate.argtypes = (c_void_p,)
def triangulate(m):
    lib.shortest_edge_triangulate(m.obj)

lib.GLManifoldViewer_new.restype = c_void_p
lib.GLManifoldViewer_delete.argtypes = (c_void_p,)
lib.GLManifoldViewer_display.argtypes = (c_void_p,c_void_p,c_char,c_bool, POINTER(c_float*3), POINTER(c_double),c_bool,c_bool)
lib.GLManifoldViewer_get_annotation_points.restype = c_size_t
lib.GLManifoldViewer_get_annotation_points.argtypes = (c_void_p, POINTER(POINTER(c_double)))
lib.GLManifoldViewer_event_loop.argtypes = (c_bool,)
class GLManifoldViewer:
    def __init__(self):
        self.obj = lib.GLManifoldViewer_new()
    def __del__(self):
        lib.GLManifoldViewer_delete(self.obj)
    def display(self, m, mode='w', smooth=True, bg_col=[0.3,0.3,0.3], data=None, reset_view=False, once=False):
        data_ct = np.array(data,dtype=c_double).ctypes
        data_a = data_ct.data_as(POINTER(c_double))
        bg_col_ct = np.array(bg_col,dtype=c_float).ctypes
        bg_col_a = bg_col_ct.data_as(POINTER(c_float*3))
        lib.GLManifoldViewer_display(self.obj,m.obj,c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
    def annotation_points(self):
        pos = POINTER(c_double)()
        n = lib.GLManifoldViewer_get_annotation_points(self.obj, byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
    @staticmethod
    def event_loop():
        lib.GLManifoldViewer_event_loop(False)

    # MeshDistance* MeshDistance_new(HMesh::Manifold* m);

lib.MeshDistance_new.restype = c_void_p
lib.MeshDistance_new.argtypes = (c_void_p,)
lib.MeshDistance_signed_distance.restype = c_float
lib.MeshDistance_signed_distance.argtypes = (c_void_p,POINTER(c_float*3),c_float)
lib.MeshDistance_delete.argtypes = (c_void_p,)
class MeshDistance:
    def __init__(self,m):
        self.obj = lib.MeshDistance_new(m.obj)
    def __del__(self):
        lib.MeshDistance_delete(self.obj)
    def signed_distance(self,p,upper=1e30):
        p_ct = np.array(p,dtype=c_float).ctypes.data_as(POINTER(c_float*3))
        return lib.MeshDistance_signed_distance(self.obj,p_ct,upper)
