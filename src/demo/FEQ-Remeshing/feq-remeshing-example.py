#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl, spatial
from glob import glob
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 
graph_dir = '../../../data/Graphs/'
# viewer = gl.Viewer()

import glob

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + 'arma*.obj')
# print(obj_files)
obj_files.sort()
for o_file in obj_files:    
    base_name = path.basename(o_file).split('.')[0]
    print("Remeshing " + o_file, flush=True)
    ref_mesh = hmesh.load(o_file)
    # viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])

    for mode in ['ps']:
        s = graph.load(graph_dir+base_name+'.graph')
        if 'p' in mode:
            graph.prune(s)
        if 's' in mode:
            graph.smooth(s, alpha=0.5)

        avg_edge_len = s.average_edge_length()
        print(avg_edge_len)
        print('Building FEQ', flush=True)
        m_skel = hmesh.skeleton_to_feq(s, node_radii=0.5*avg_edge_len, symmetrize=True)
        hmesh.save(base_name + "-" + mode + "-skel.obj", m_skel)
        print('Fitting to reference mesh')
        fit_mesh = hmesh.Manifold(m_skel)
        hmesh.cc_subdivide(fit_mesh)
        hmesh.fit_mesh_mmh(fit_mesh, ref_mesh, iterations=150)
        # hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, iter=15)
        # hmesh.cc_subdivide(fit_mesh)
        # hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, iter=25)

        # for iter in range(5):
        #     for _ in range(10):
        #         fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, iter=1)
        #         hmesh.cc_smooth(fit_mesh)
        #     if iter < 2:
        #         hmesh.cc_subdivide(fit_mesh)

        # viewer.display(fit_mesh, bg_col=[1,1,1], reset_view=True)
        out_file = base_name + "-" + mode + "-out.obj"
        hmesh.save(out_file, fit_mesh)

