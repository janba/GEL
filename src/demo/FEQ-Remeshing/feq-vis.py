#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl, spatial
from glob import glob
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 
graph_dir = '../../../data/Graphs/'
viewer = gl.Viewer()

import glob

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + '*.obj')
obj_files.sort()
for o_file in obj_files:    
    base_name = path.basename(o_file).split('.')[0]
    ref_mesh = hmesh.load(o_file)
    s = graph.load(graph_dir+base_name+'.graph')
    if base_name != "warrior":
        graph.prune(s)
        graph.smooth(s, alpha=0.5)
        out_mesh = hmesh.load(base_name + "-ps-out.obj")
    else:
        out_mesh = hmesh.load(base_name + "--out.obj")

    viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])
    viewer.display(ref_mesh, g=s, mode='x', bg_col=[1,1,1])
    viewer.display(out_mesh, reset_view=False, bg_col=[1,1,1])
