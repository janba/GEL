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
    print("Remeshing " + o_file)
    ref_mesh = hmesh.load(o_file)
    # viewer.display(ref_mesh, reset_view=True, bg_col=[1,1,1])

    for mode in ['', 'ps']:
# for lap_wt in [0.2, 0.25, 0.35, 0.5, 1]:
        # mode = ''
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
        # for _ in range(1):
        #     hmesh.non_rigid_registration(fit_mesh, ref_mesh)
        #     # hmesh.regularize_quads(fit_mesh, w=0.75, shrink=0.3, iter=3)
        #     print(".")

        # ref_mesh_orig = hmesh.Manifold(ref_mesh)
        fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh)
        # viewer.display(fit_mesh, bg_col=[1,1,1])
        # ref_mesh = hmesh.fit_mesh_to_ref(ref_mesh, fit_mesh, dist_wt=0.5, lap_wt=1, iter=1)
        # viewer.display(ref_mesh, bg_col=[1,1,1])

        # ref_tree = spatial.I3DTree()
        # ref_pos = ref_mesh.positions()
        # for v in ref_mesh.vertices():
        #     ref_tree.insert(ref_pos[v],v)
        # ref_tree.build()

        # fit_pos = fit_mesh.positions()
        # ref_pos_orig = ref_mesh_orig.positions()
        # for v in fit_mesh.vertices():
        #     p = fit_pos[v]
        #     _, ref_v = ref_tree.closest_point(p, 1e32)
        #     # I = [ref_v] + list(ref_mesh.circulate_vertex(ref_v))
        #     # R,T = least_squares_affine_transformation(ref_pos[I], ref_pos_orig[I])
        #     fit_pos[v] = (ref_pos_orig[ref_v]-ref_pos[ref_v]) + p

        viewer.display(fit_mesh, bg_col=[1,1,1], reset_view=True)
        out_file = base_name + "-" + mode + "-out.obj"
        hmesh.save(out_file, fit_mesh)

