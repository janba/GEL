#!/opt/local/bin/python
#
# This script computes the set of homology generators for a mesh surface.
# Subsequently it cuts the mesh along these curves, producing a new mesh that is
# topologically equivalent to a disk. Note that the script does not handle 
# large meshes well, so it reduces the mesh using quadric simplification.
# The script also closes holes in the mesh before computing the homology generators.
#
# Usage: python homology_generators.py <mesh_file>

from numpy import zeros, array, argmax
from numpy.linalg import norm
from scipy.spatial import KDTree
from pygel3d import graph, hmesh, gl_display as gl
from queue import Queue, PriorityQueue
from sys import argv

def graph_from_mesh_verts(m: hmesh.Manifold):
    ''' This simple function converts a mesh, m, to a graph with a set of nodes 
    that precisely corresponds to the set of vertices of m and with the same numbering.
    No edges are added to the graph. '''
    g_out = graph.Graph()
    pos = m.positions()
    for v in m.vertices():
        g_out.add_node(pos[v])
    return g_out

def mst(m: hmesh.Manifold, v0=None):
    ''' Given a mesh, m, as the argument and an optional starting vertex, v0, mst 
    returns the minimum spanning tree of g. The spanning tree is encoded as a vector
    of halfedge tags. A halfedge is tagged if it belongs to the MST and points out 
    from the root.'''
    if v0 is None:
        v0 = m.vertices()[0]
    htag = zeros(m.no_allocated_halfedges())
    pos = m.positions()
    edges = []
    for v in m.vertices():
        for h in m.circulate_vertex(v, mode='h'):
            vn = m.incident_vertex(h)
            edges += [ (norm(pos[vn]-pos[v]), v, vn, h) ]
    edges.sort()
    visited = { v0 }
    while edges:
        for _, i, j, h in edges:
            if i in visited and not j in visited:
                htag[h] = 1
                visited.add(j)
                break
        edges = [ (d, i, j, h) for d, i, j, h in edges if not j in visited ]
    return htag

def shortest_path_tree(m: hmesh.Manifold, v0=None):
    ''' Given a mesh, m, as the argument and an optional starting vertex, v0, 
    shortest_path_tree returns the tree of shartest paths from v0 to all other
    vertices in the mesh. The spanning tree is encoded as a set of halfedge tags.
    A halfedge is tagged if it belongs to the spanning tree and points out. The 
    distances are also returned as a vector with an entry per vertex.'''
    dist = array([1e300]*m.no_allocated_vertices())
    if v0 is None:
        v0 = m.vertices()[0]
    htag = zeros(m.no_allocated_halfedges())
    pos = m.positions()
    q = PriorityQueue()
    q.put((0, v0, None))
    while not q.empty():
        d, v, h = q.get()
        if d < dist[v]:
            if not h is None:
                htag[h] = 1
            dist[v] = d
            for hn in m.circulate_vertex(v, mode='h'):
                vn = m.incident_vertex(hn)
                dn = d + norm(pos[vn]-pos[v])
                if dn < dist[vn]:
                    q.put((dn, vn, hn))
    return htag, dist


def find_instigator_edges(m: hmesh.Manifold, htag, f0=None):
    '''This function is passed a mesh m and a set of halfedge tags that define the
    MST. Instigator edges are edges that do not belong to the MST but which form loops that
    cut the surface without dividing it into two disjoint parts. We find these instigators
    by creating a spanning tree on the dual graph of the mesh. If two adjacent faces (neither
    being the other's parent) are both incident to an edge that is not in the MST, then 
    the two faces are connected by an instigator edge. The instigator edges are returned 
    as a set of halfedges. '''
    F = m.faces()
    if f0 is None:
        f0 = F[0]
    touched = zeros(m.no_allocated_faces())
    touched[f0] = 1
    Q = Queue()
    Q.put((f0,f0))
    instigators = set()
    while not Q.empty():
        f, f_parent = Q.get()
        for h in m.circulate_face(f, mode='h'):
            ho = m.opposite_halfedge(h)
            fo = m.incident_face(ho)
            if fo != f_parent:
                if not (htag[h]==1 or htag[ho]==1):
                    if touched[fo]==1:
                        instigators.add(min(h,ho))
                    else:
                        Q.put((fo,f))
                        touched[fo] = 1
    return instigators

def trace_back(m: hmesh.Manifold, htag, v):
    '''Given a mesh, m, a spanning tree, htag, in the form of a set of halfedge tags
    and a starting vertex, v, this code traces back from v to the root of the
    spanning tree, producing a path which is returned as a list of halfedges.'''
    l = []
    while v:
        found = False
        for h in m.circulate_vertex(v, mode='h'):
            ho = m.opposite_halfedge(h)
            if htag[ho]==1:
                l.append(h)
                v =  m.incident_vertex(h)
                found = True
                break
        if not found:
            v = None
    return l

def form_loop(m: hmesh.Manifold, htag, i):
    '''Given a mesh, m, and a spanning tree, htag, in the form of a set of halfedge tags, htag, 
    and a mesh edge i (instigator), we trace back from both end points to the spanning 
    tree root, and combine the two curves to form a loop. The loop is returned as a list of 
    successive halfedges.'''
    loop0 = trace_back(m, htag, m.incident_vertex(i))
    loop1 = trace_back(m, htag, m.incident_vertex(m.opposite_halfedge(i)))
    loop  = [ m.opposite_halfedge(h) for h in loop0 if not h in loop1 ]
    loop += [ i ]
    loop += [ h for h in loop1 if not h in loop0 ]
    return loop

def loop_to_graph(m: hmesh.Manifold, loop):
    '''Given a mesh and a loop, this function creates a graph with the same
    number of nodes as the number of vertices in the loop. The edges of the graph
    are the edges of the loop. The graph is returned.'''
    g_loop = graph_from_mesh_verts(m)
    for h in loop:
        v0 = m.incident_vertex(h)
        v1 = m.incident_vertex(m.opposite_halfedge(h))
        g_loop.connect_nodes(v0, v1)
    for n in g_loop.nodes():
        if len(g_loop.neighbors(n)) == 0:
            g_loop.remove_node(n)
    return g_loop

def cut_mesh(m: hmesh.Manifold, loops):
    ''' Given a mesh, m, and a list of loops, this function
    cuts the mesh along the edges defined by loops. The cut is done by
    creating a new mesh, m_out, and copying the vertices and faces of m
    to m_out. The edges of m_out are then stitched together except where the
    edges are cut. The function returns the new mesh, m_out. '''
    cut_tags = zeros(m.no_allocated_halfedges())
    for loop in loops:
        for h in loop:
            ho = m.opposite_halfedge(h) 
            cut_tags[h] = 1
            cut_tags[ho] = 1
    m_out = hmesh.Manifold()
    h_map = {}
    pos = m.positions()
    for f in m.faces():
        f_out = m_out.add_face(pos[m.circulate_face(f)])
        h_m = m.circulate_face(f, mode='h')
        h_out = [ m_out.opposite_halfedge(h) for h in m_out.circulate_face(f_out, mode='h') ]
        for i, h in enumerate(h_m):
            h_map[h] = h_out[i]

    for h in m.halfedges():
        if cut_tags[h]==0:
            ho = m.opposite_halfedge(h)
            if h < ho:
                h0 = h_map[h]
                h1 = h_map[ho]
                h0_in_use = m_out.halfedge_in_use(h0)
                h1_in_use = m_out.halfedge_in_use(h1)
                if h0_in_use and h1_in_use:
                    m_out.stitch_boundary_edges(h0, h1)
                else:
                    print("Should not happen: one or more halfedges unexpectedly not in use")
    m_out.cleanup()
    return m_out


if __name__ == "__main__":
    m = hmesh.load(argv[1])
    if m is None:
        print("Failed to load mesh.")
        exit(1)
    # hmesh.quadric_simplify(m,0.25)
    hmesh.close_holes(m, 100000)
    m.cleanup()

    # We now compute the genus directly from the number of vertics, edges
    # and faces using the Euler-Poincare formula:
    #    V - E + F = 2 (s-g) - b
    # where V is the number of vertices, E number of edges, and F is the
    # number of faces:
    V = len(m.vertices())
    E = len(m.halfedges())//2
    F = len(m.faces())

    # s is the number of shells, b is the number of boundary
    # curves, and g is the genus. We assume b is zero and s is 1. This is
    # tantamount to assuming that the surfaces is a single closed mesh.
    # Hence, starting from:
    #    V - E + F = 2 - 2 g
    # we can solve for g, obtaining:
    g = (2 - (V - E + F))//2
    print("genus: ", g)

    # Next, we use Crane's method (https://www.cs.cmu.edu/~kmcrane/Projects/LoopsOnSurfaces/)
    # to find a set of curves that cut the surface into a disk. We can think of these curves
    # as a basis for the homology group H1 or the fundamental group of the surface.
    # The number of such curves is two times the genus. Hence, this gives us an independent 
    # way of finding the genus in addition to a way of cutting the surface.
    # Display the cut curves
    viewer = gl.Viewer()
    viewer.display(m)

    # We need to find a starting point for the MST. We can let the user do
    # this by clicking on a vertex. We use a KDTree to find the closest vertex
    # to the clicked point. Alternatively, we can just start from the first
    # vertex in the mesh.
    tree = KDTree(m.positions())
    try:
        a_pts = viewer.annotation_points()
        v0 = tree.query(a_pts[0])[1]
    except:
        # If we can't find the vertex, we just use the first one.
        v0 = m.vertices()[0]

    htag, dist = shortest_path_tree(m, v0)

    f0 = m.circulate_vertex(argmax(dist), mode='f')[0]
    inst = find_instigator_edges(m, htag)
    print("genus from H1 basis: ", len(inst)//2)

    loops = []
    for i in inst:
        loop = form_loop(m,htag,i)
        loops.append(loop)
        viewer.display(m, loop_to_graph(m, loop), mode="g", smooth=False)

    m_cut = cut_mesh(m, loops)
    hmesh.save(f"{argv[1].split('.')[0]}_cut.obj", m_cut)
    hmesh.laplacian_smooth(m_cut, 0.5, 1)
    viewer.display(m_cut, smooth=True, mode='g')
