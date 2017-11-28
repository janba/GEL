from ctypes import *
#from numpy.ctypeslib import as_array
import numpy as np
import random
from copy import copy
from numpy.linalg import norm
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
    def __init__(self,orig=None):
        if orig == None:
            self.obj = lib.Manifold_new()
        else:
            self.obj = lib.Manifold_copy(orig.obj)
    def __del__(self):
        lib.Manifold_delete(self.obj)
    def add_face(self,pts):
        lib.Manifold_add_face(self.obj, len(pts), np.array(pts))
    def positions(self):
        pos = POINTER(c_double)()
        n = lib.Manifold_positions(self.obj, byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
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

lib.GLManifoldViewer_new.restype = c_void_p
lib.GLManifoldViewer_delete.argtypes = (c_void_p,)
lib.GLManifoldViewer_display.argtypes = (c_void_p,c_void_p,c_char,c_bool, POINTER(c_float*3), POINTER(c_double),c_bool)
class GLManifoldViewer:
    def __init__(self):
        self.obj = lib.GLManifoldViewer_new()
    def __del__(self):
        lib.GLManifoldViewer_delete(self.obj)
    def display(self, m, mode='w', smooth=True, bg_col=[0.3,0.3,0.3],data=None, once=False):
        data_ct = np.array(data,dtype=c_double).ctypes
        data_a = data_ct.data_as(POINTER(c_double))
        bg_col_ct = np.array(bg_col,dtype=c_float).ctypes
        bg_col_a = bg_col_ct.data_as(POINTER(c_float*3))
        lib.GLManifoldViewer_display(self.obj,m.obj,c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,once)


def test_I3DTree():
    t = I3DTree()
    t.insert([1,1,1],1)
    t.insert([2,2,2],2)
    t.insert([2,1,2],3)
    t.insert([1,2,1],4)
    t.insert([1,1,2],5)
    t.insert([2,1,1],6)
    t.build()

    print("Found point: " , t.closest_point([1.2,1.9,0.9],0.4))

    (K,V) = t.in_sphere([1.2,1.9,0.9],1.5)

    for k,v in zip(K,V):
        print(k, v)

def laplacian_smooth(m, max_iter=1):
    pos = m.positions()
    for iter in range(0,max_iter):
        old_pos = pos[:,:]
        for v0 in m.vertices():
            p0 = old_pos[v0]
            avg_p = [.0,.0,.0]
            w_sum = 0
            for vi in m.circulate_vertex(v0,'v'):
                avg_p += old_pos[vi]
                w_sum += 1.0
            pos[v0] = avg_p / w_sum

def dual(m):
    pos = m.positions()
    m2 = Manifold()
    for v in m.vertices():
        pts = []
        for f in m.circulate_vertex(v,'f'):
            avg_p = [.0,.0,.0]
            w_sum = 0
            for vf in m.circulate_face(f, 'v'):
                avg_p += pos[vf]
                w_sum += 1
            pts += [avg_p/w_sum]
        m2.add_face(pts)
    stitch(m2)
    return m2

def test_GLManifoldViewer():
    m = obj_load("/Users/janba/Teaching/02580/pyigl-experiments/bunnygtest.obj")
    m2 = dual(m)
    pos = m.positions()
    viewer = GLManifoldViewer()
    vfield = pos[:,:]
    viewer.display(m,'l',smooth=True,data=vfield,bg_col=[0.9,0.3,0.3])

    sfield = pos[:,1]
    viewer.display(m,'s',smooth=True,data=sfield)
    viewer.display(m,'g',smooth=True)
    viewer.display(m,'n',smooth=False)

    viewer.display(m,'i',smooth=True,data=sfield)
    laplacian_smooth(m,10)
    viewer.display(m,'i',smooth=True,data=sfield)

    viewer.display(m2,'w',smooth=True,data=sfield)

def average_edge_length(m):
    hedges = m.halfedges()
    lsum = 0
    for h in hedges:
        lsum += m.edge_length(h)
    return lsum / len(hedges)

def delaunay_edge(m,h):
    if m.is_halfedge_at_boundary(h):
        return True
    pos = m.positions()
    mat = np.ones((4,4))
    i = 0
    for v in m.circulate_face(m.incident_face(h)):
        mat[i,:] = [pos[v,0],pos[v,1],pos[v,0]**2 + pos[v,1]**2,1]
        i += 1
    v = m.incident_vertex(m.next_halfedge(m.opposite_halfedge(h)))
    mat[3,:] = [pos[v,0],pos[v,1],pos[v,0]**2 + pos[v,1]**2,1]
    return True if np.linalg.det(mat)<1e-300 else False

def fractal_terrain():
    m = Manifold()
    m.add_face([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    m.add_face([[1.,1.,0.],[0.,1.,0.],[1.,0.,0.]])
    stitch(m)
    for iter in range(0,5):
        for f in m.faces():
            v_new = m.split_face_by_vertex(f)
            m.positions()[v_new,2] += 0.3 * 0.5**iter * random.random()
        for h in m.halfedges():
            if not delaunay_edge(m,h):
                m.flip_edge(h)

    # Test remove_edge
    avg_len = average_edge_length(m)
    for h in m.halfedges():
        if m.halfedge_in_use(h):
            if m.edge_length(h)>2.0 * avg_len:
                m.remove_edge(h)

    Pmin,Pmax = bbox(m)
    print("Bbox corners ", Pmin[:], Pmax[:])
    C,r = bsphere(m)
    print("Center =", C[:]," rad = ", r)
    return m

def test_Manifold_class():

    m = fractal_terrain()
    m_backup = Manifold(m)
    viewer = GLManifoldViewer()
    viewer.display(m,'w',smooth=False)

    # Test Manifold_merge_faces
    for h in m.halfedges():
        if m.halfedge_in_use(h) and random.random() > 0.9:
            m.merge_faces(h)
    viewer.display(m,'w',smooth=False)

    m = Manifold(m_backup)

    # Test remove_vertex
    for v in m.vertices():
        if random.random() > 0.95:
            m.remove_vertex(v)
    viewer.display(m,'w',smooth=False)
    # do cleanup validate that the shape of positions array shrinks
    print(m.positions().shape)
    m.cleanup()
    print(m.positions().shape)
    # Test remove_face
    for f in m.faces():
        if random.random() > 0.95:
            m.remove_face(f)
    viewer.display(m,'w',smooth=False)
    # Test validity
    print("valid : ", valid(m))
    # Test closed
    print("closed : ", closed(m))
    for h in m.halfedges():
        if not m.face_in_use(m.incident_face(h)):
            m.close_hole(h)
    print("closed : ", closed(m))
    viewer.display(m,'w',smooth=False)

    m = Manifold(m_backup)

    # Test Manifold collapse_edge
    for h in m.halfedges():
        if m.halfedge_in_use(h) and random.random() > 0.9:
            m.collapse_edge(h)
    viewer.display(m,'w',smooth=False)

def make_tet():
    verts = [[ 0.0000000e+0, 2.00000000, 0.0000000e+0],
     [0.0000000e+0 ,-0.66666000 ,1.88561800],
     [-1.63299400, -0.66666600, -0.94281000],
     [1.63299400, -0.66666600, -0.94281000]]
    m = Manifold()
    m.add_face([verts[0],verts[2],verts[1]])
    m.add_face([verts[0],verts[3],verts[2]])
    m.add_face([verts[1],verts[3],verts[0]])
    m.add_face([verts[2],verts[3],verts[1]])
    stitch(m)
    m.cleanup()
    return m

def make_grid(Ni,Nj,scale=1.0):
    m = Manifold()
    for i in range(0,Ni):
        for j in range(0,Nj):
            x = float(i)*scale
            y = float(j)*scale
            m.add_face([[x,y,.0],[x+scale,y,.0],[x+scale,y+scale,.0],[x,y+scale,.0]])
    stitch(m)
    m.cleanup()
    return m

def test_dynamics():
    mass = 1.0
    k = 1.0
    dt = 0.5
    m = make_grid(30,20)
    V = m.vertices()
    pos = m.positions()
    H = m.halfedges()
    NH = len(H)
    NV = len(V)

    F_ext = np.zeros((NV,3))

    load_vertices = []
    support_vertices = []
    for v in V:
        if norm(pos[v]-[3.0,0.0,0.0])<1e-10:
            load_vertices += [v]
        if pos[v,0]<1e-10:
            support_vertices += [v]


    rest_length = np.zeros(NH)
    for h in m.halfedges():
        rest_length[h] = m.edge_length(h)
    viewer = GLManifoldViewer()
    velocity = np.zeros((NV,3))

    for v in load_vertices:
        F_ext[v] = [0,-0.005,0]


    for iter in range(0,10000):
        for v in V:
            F = np.copy(F_ext[v])
            # if v in load_vertices:
            #     F = [.0,-1,.0]#F_ext[v]
            for h in m.circulate_vertex(v,mode='h'):
                vn = m.incident_vertex(h)
                l = m.edge_length(h)
                rl = rest_length[h]
                d = (pos[vn]-pos[v])/l
                F += k*(l-rl)*d
            velocity[v] = 0.9 * velocity[v] + dt * F/mass
        for v in V:
            if v not in support_vertices:
                pos[v] += dt * velocity[v]
#            print(pos[0,:], pos[1,:], pos[2,:], pos[3,:])
        viewer.display(m,'w',once=True)


class SphereVolume:
    def __init__(self,_center,_radius):
        self.center = _center
        self.radius = _radius
    def eval(self,p):
        return norm(p-self.center) - self.radius
    def grad(self,p):
        v = p - self.center
        return v / norm(v)

class XForm:
    def __init__(self,Plow,Phigh,dims):
        self.dims = np.array(dims)
        self.Plow = np.array(Plow)
        self.Phigh = np.array(Phigh)
        self.scaling = (self.Phigh - self.Plow) / self.dims
        self.inv_scaling = self.dims / (self.Phigh - self.Plow)
    def map(self,pi):
        return self.scaling * np.array(pi) + self.Plow
    def inv_map(self,p):
        return floor(np.array(p-Plow) * self.inv_scaling)


cube_faces = np.array([[[0,-0.5,-0.5],[0,0.5,-0.5],[0,0.5,0.5],[0,-0.5,0.5]],
    [[0, 0.5,-0.5],[0,-0.5,-0.5],[0,-0.5,0.5],[0,0.5,0.5]],
    [[ 0.5,0, -0.5],[-0.5,0, -0.5],[-0.5,0, 0.5],[0.5,0, 0.5]],
    [[-0.5,0, -0.5],[0.5,0, -0.5],[0.5,0, 0.5],[-0.5,0, 0.5]],
    [[-0.5,-0.5,0],[0.5,-0.5,0],[0.5,0.5,0],[-0.5,0.5,0]],
    [[ 0.5,-0.5,0],[-0.5,-0.5,0],[-0.5,0.5,0],[0.5,0.5,0]]])

neighbours = np.array([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.],[0.,0.,1.],[0.,0.,-1.]])

def polygonize_volume(vol,xform,tau):
    m = Manifold()
    for pi in np.ndindex(tuple(xform.dims)):
        if vol.eval(xform.map(pi)) < tau:
            for l in range(0,6):
                nbr = neighbours[l]
                p_nbr = xform.map(pi+nbr)
                if vol.eval(p_nbr) >= tau:
                    m.add_face(cube_faces[l]+pi+0.5*nbr)
    stitch(m)
    m.cleanup()
    pos = m.positions()
    for v in m.vertices():
        p = xform.map(pos[v])
        for _ in range(0,3):
            g = vol.grad(p)
            p -= vol.eval(p) * g / np.dot(g,g)
        pos[v] = p

    return m

def test_volume_polygonize():
    vol = SphereVolume([0.,0.5,0.],1.4)
    xform = XForm([-1.5,-1.5,-1.5],[1.5,1.5,1.5],[16,16,16])
    m = polygonize_volume(vol, xform, 0)
    viewer = GLManifoldViewer()
    viewer.display(m,'w')

#test_volume_polygonize()
#test_dynamics()
#test_Manifold_class()
test_GLManifoldViewer()
#test_I3DTree()
