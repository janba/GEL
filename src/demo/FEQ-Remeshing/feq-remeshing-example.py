#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl
from glob import glob
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 

viewer = gl.Viewer()

import glob

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + '*.obj')
for o_file in obj_files:
    print("Remeshing " + o_file)
    ref_mesh = hmesh.load(o_file)
    viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])

    print('Building skeleton')
    g = graph.from_mesh(ref_mesh)
    s = graph.MSLS_skeleton(g, grow_thresh=512)
    for _ in range(3):
        graph.smooth(s)
    
    print('Building FEQ')
    m_skel = hmesh.skeleton_to_feq(s, symmetrize=True)
    viewer.display(m_skel, bg_col=[1,1,1])

    print('Refining FEQ')
    hmesh.cc_split(m_skel)
    for _ in range(5):
        hmesh.cc_smooth(m_skel)
    viewer.display(m_skel, bg_col=[1,1,1])

    print('Fitting to reference mesh')
    fit_mesh = hmesh.Manifold(m_skel)
    fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, dist_wt=1, lap_wt=1)
    viewer.display(fit_mesh, bg_col=[1,1,1])
    out_file = path.basename(o_file) + "-out.obj"
    hmesh.save(out_file, fit_mesh)

