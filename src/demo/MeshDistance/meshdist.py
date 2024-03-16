from pygel3d import  hmesh
from numpy import array, vstack

m = hmesh.load("../../../data/as.obj")
a = array(list(hmesh.bbox(m)) + [array((0.0,20.0,0.0))])
print(a)
D = hmesh.MeshDistance(m)

print("d(a[0])", D.signed_distance(a[0]))
print("d(a[1])", D.signed_distance(a[1]))
b = [array([a[0,0], a[1,0]]), array([a[0,1], a[1,1]]), array([a[0,2], a[1,2]])]
c = vstack(b).T
d = c.flatten()
print("a.flags", a.flags)
print("c.flags", c.flags)
print("d(a)", D.signed_distance(a))
print("d(c)", D.signed_distance(c))
print("d(d)", D.signed_distance(d))
print("I(a)", D.ray_inside_test(a))
print("I(c)", D.ray_inside_test(c))
print("I(d)", D.ray_inside_test(d))

for v in a:
    print(D.intersect(v, array((0.0, 1.0, 0.0))))