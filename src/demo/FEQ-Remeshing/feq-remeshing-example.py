#!/opt/local/bin/python
from pygel3d import hmesh, graph, gl_display as gl, spatial
from glob import glob
from numpy import zeros, exp, array
from os import path
mesh_dir = '../../../data/ReferenceMeshes/' 
graph_dir = '../../../data/Graphs/'
# viewer = gl.Viewer()

import glob

def edge_len_ratio(m: hmesh.Manifold, v):
    l = []
    for h in m.circulate_vertex(v,'h'):
        l.append(m.edge_length(h))
    return 0.2# if len(l)==4 else 0.2

mesh_dir = '../../../data/ReferenceMeshes/'
obj_files = glob.glob(mesh_dir + '*.obj')
obj_files.sort()
for o_file in obj_files:    
    base_name = path.basename(o_file).split('.')[0]
    print("Remeshing " + o_file, flush=True)
    ref_mesh = hmesh.load(o_file)
    for mode in ['ps']:
        s = graph.load(graph_dir+base_name+'.graph')
        if 'p' in mode:
            graph.prune(s)
        if 's' in mode:
            graph.smooth(s, alpha=0.5)

        radii = zeros(len(s.nodes()))
        pos = s.positions()
        avg_edge_len = s.average_edge_length()
        for n in s.nodes():
            radii[n] = 0.5*avg_edge_len
        print(avg_edge_len)
        print('Building FEQ', flush=True)
        m_skel = hmesh.skeleton_to_feq(s, node_radii=radii, symmetrize=True)
        hmesh.save(base_name + "-" + mode + "-skel.obj", m_skel)
        print('Fitting to reference mesh')
        fit_mesh = hmesh.Manifold(m_skel)
        mesh_dist = hmesh.MeshDistance(ref_mesh)
        print('Stable marriage registration', end='')
        for i in range(2):
            print('.',end='',flush=True)
            hmesh.cc_subdivide(fit_mesh)
            hmesh.stable_marriage_registration(fit_mesh, ref_mesh)
        print('')
        print('Inflating and adjusting', end='')
        for i in range(12):
            print('.',end='',flush=True)
            hmesh.inflate_mesh(fit_mesh, mesh_dist=mesh_dist)
            hmesh.volume_preserving_cc_smooth(fit_mesh, iter=5)
            hmesh.kill_degenerate_face_loops(fit_mesh, thresh=0.25)
            hmesh.volume_preserving_cc_smooth(fit_mesh, iter=5)
        out_file = base_name + "-" + mode + "-out.obj"
        hmesh.save(out_file, fit_mesh)

