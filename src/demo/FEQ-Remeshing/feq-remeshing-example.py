#!/opt/local/bin/python
from this import d
from pygel3d import hmesh, graph, gl_display as gl
from os import getcwd

graphs = [
'hand.graph',
'armadillo_symmetric.graph',
'bunny.graph',
'feline.graph',
'fertility.graph',
'warrior.graph']

objs = [
'usai_hand_tri.obj',
'armadillo.obj',
'bunny.obj',
'feline.obj',
'fertility_tri.obj',
'warrior.obj'
]

iters = [(250,0.9,1.0), (50,0.5,1.0), (50,0.5,1.0), (50,0.5,1.0), (50,0.5,1.0), (50,0.5,1.0)]

mesh_dir = '../../../data/ReferenceMeshes/' 
skel_dir = '../../../data/Graphs/'

viewer = gl.Viewer()

for g_file, o_file, params in zip(graphs, objs, iters):
    iter, dist_wt, lap_wt = params
    print("Remeshing " + o_file)

    print('Building FEQ')
    s = graph.load(skel_dir + g_file)
    m_skel = hmesh.skeleton_to_feq(s)#, [5.0]*len(s.nodes()))
    hmesh.cc_split(m_skel)
    hmesh.cc_smooth(m_skel)

    print('Fitting to reference mesh')
    ref_mesh = hmesh.load(mesh_dir + o_file)
    fit_mesh = hmesh.Manifold(m_skel)
    fit_mesh = hmesh.fit_mesh_to_ref(fit_mesh, ref_mesh, local_iter=iter, dist_wt=dist_wt, lap_wt=lap_wt)

    print("Displaying. HIT ESC IN GRAPHICS WINDOW TO PROCEED...")
    viewer.display(fit_mesh, reset_view=True)
    hmesh.save(o_file + "-out.obj", fit_mesh)
    # viewer.display(m_skel, reset_view=True)

