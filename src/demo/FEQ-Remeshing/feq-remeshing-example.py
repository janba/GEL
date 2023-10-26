#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl
from glob import glob
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 
graph_dir = '../../../data/Graphs/'
# viewer = gl.Viewer()

import glob

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + 'armadillo.obj')
obj_files.sort()
for o_file in obj_files:    
    base_name = path.basename(o_file).split('.')[0]
    print("Remeshing " + o_file)
    ref_mesh = hmesh.load(o_file)
    # viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])

#    for mode in ['', 'ps']:
for lap_wt in [0.1, 0.25, 0.5, 0.75, 1]:
        mode = ''
        s = graph.load(graph_dir+base_name+'.graph')
        if 'p' in mode:
            graph.prune(s)
        if 's' in mode:
            graph.smooth(s, alpha=0.5)

        print('Building FEQ')
        m_skel = hmesh.skeleton_to_feq(s, node_radii=0.0, symmetrize=True)
        # viewer.display(m_skel, bg_col=[1,1,1])

        print('Fitting to reference mesh')
        fit_mesh = hmesh.Manifold(m_skel)
        fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, dist_wt=1, lap_wt=lap_wt, iter=10)
        # hmesh.cc_split(fit_mesh)
        # fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, dist_wt=1, lap_wt=1, iter=10)
        # viewer.display(fit_mesh, bg_col=[1,1,1])
        out_file = base_name + "-" + str(lap_wt) + "-out.obj"
        hmesh.save(out_file, fit_mesh)

