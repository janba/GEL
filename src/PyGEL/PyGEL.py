from ctypes import *
lib = cdll.LoadLibrary("libPyGEL.dylib")

Vec3d = c_double * 3

lib.IntVector_new.restype = c_void_p
lib.IntVector_get.argtypes = (c_void_p, c_int)
lib.IntVector_delete.argtypes = (c_void_p,)
class IntVector:
    def __init__(self,data,s):
        self.obj = data
        self.size = s
    def __del__(self):
        lib.IntVector_delete(self.obj)
    def __iter__(self):
        for i in range(0,self.size):
            yield lib.IntVector_get(self.obj, i)

lib.Vec3dVector_new.restype = c_void_p
lib.Vec3dVector_get.argtypes = (c_void_p, c_int)
lib.Vec3dVector_get.restype = POINTER(c_double)
lib.Vec3dVector_delete.argtypes = (c_void_p,)
class Vec3dVector:
    def __init__(self,data,s):
        self.obj = data
        self.size = s
    def __del__(self):
        lib.Vec3dVector_delete(self.obj)
    def __iter__(self):
        for i in range(0,self.size):
            data = lib.Vec3dVector_get(self.obj, i)
            yield [data[0], data[1], data[2]]


lib.I3DTree_new.restype = c_void_p
lib.I3DTree_delete.argtypes = (c_void_p,)
lib.I3DTree_insert.argtypes = (c_void_p, c_double, c_double, c_double, c_int)
lib.I3DTree_build.argtypes = (c_void_p,)
lib.I3DTree_closest_point.argtypes = (c_void_p, c_double, c_double, c_double, c_double, POINTER(Vec3d), POINTER(c_int))
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
        key = Vec3d()
        val = c_int()
        n = lib.I3DTree_closest_point(self.obj, p[0],p[1],p[2],r,pointer(key),pointer(val))
        if n==1:
            return ([key[0],key[1],key[2]],val.value)
        return None
    def in_sphere(self, p, r):
        keys = lib.Vec3dVector_new(0)
        vals = lib.IntVector_new(0)
        n = lib.I3DTree_in_sphere(self.obj, p[0],p[1],p[2],r,keys,vals)
        return (Vec3dVector(keys,n),IntVector(vals,n))

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

test_I3DTree()
