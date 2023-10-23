#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl
from glob import glob
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 
graph_dir = '../../../data/Graphs/'
viewer = gl.Viewer()

import glob

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + '*.obj')
obj_files.sort()
obj_files.reverse()
# for o_file in [mesh_dir+'warrior.obj',mesh_dir+'torus.obj']:
# for o_file in [mesh_dir+'torus.obj']:
for o_file in obj_files:    
    #[mesh_dir+'torus.obj']: # obj_files: ##[mesh_dir+'wolf.obj', mesh_dir+'bunny.obj']: #    
    base_name = path.basename(o_file).split('.')[0]
    print("Remeshing " + o_file)
    ref_mesh = hmesh.load(o_file)
    viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])

    for mode in ['']:#, 'ps']:
        s = graph.load(graph_dir+base_name+'.graph')
        if 'p' in mode:
            graph.prune(s)
        if 's' in mode:
            graph.smooth(s, alpha=0.5)

        print('Building FEQ')
        m_skel = hmesh.skeleton_to_feq(s, node_radii=0.0, symmetrize=True)
        viewer.display(m_skel, bg_col=[1,1,1])

        print('Fitting to reference mesh')
        fit_mesh = hmesh.Manifold(m_skel)
        fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, dist_wt=1, lap_wt=0.1, local_iter=10)
        viewer.display(fit_mesh, bg_col=[1,1,1])
        out_file = base_name + "-out.obj"
        hmesh.save(out_file, fit_mesh)

