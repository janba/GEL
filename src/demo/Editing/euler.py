''' This script removes faces from a mesh until no more faces can be removed.
The mesh is loaded from a file given as the first argument to the script.
Faces are removed if they have two boundary halfedges or a single boundary. 
The point of the script is that it illustrates the proof that V-E+F=2 for 
genus 0 meshes by turning them into topological discs and then reducing them
to a single face. The script needs to be called with a mesh that has disc topology
and is triangulated.'''
from pygel3d import hmesh, gl_display as gl
from sys import argv as args
from time import sleep

m = hmesh.load(args[1] if len(args)>1 else '../../../data/cube-1side.obj')

def removable(m, f):
    '''Returns True if the face f is removable, False otherwise. 
    A face is removable if it has two boundary halfedges or a single
    boundary edge and two boundary vertices.'''
    cnt = 0
    for h in m.circulate_face(f, mode='h'):
        if m.is_halfedge_at_boundary(h):
            cnt += 1
    if cnt == 2:
        return True
    if cnt == 0:
        return False
    cnt = 0
    for v in m.circulate_face(f, mode='v'):
        if m.is_vertex_at_boundary(v):
            cnt += 1
    return cnt == 2

v = gl.Viewer()
v.display(m, bg_col=[1,1,1])
did_work = True
while did_work:
    did_work = False
    for f in m.faces():
        if removable(m, f):
            m.remove_face(f)
            did_work = True
            break
    v.display(m, once=True, bg_col=[1,1,1])
    sleep(.5)
v.display(m, bg_col=[1,1,1])
del v