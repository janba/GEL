""" This module provides a Graph class and functionality for skeletonization using graphs. """
import ctypes as ct
import numpy as np
from pygel3d import hmesh, lib_py_gel, IntVector

class Graph:
    """ This class is for representing graphs embedded in 3D. The class does not in
    itself come with many features: it contains methods for creating, accessing, and
    housekeeping. When vertices are used as parameters in the functions below, we usually
    use the parameter name n (for node). n is simply an index (i.e. an integer) that
    refers to a node (aka vertex)."""
    def __init__(self,orig=None):
        if orig == None:
            self.obj = lib_py_gel.Graph_new()
        else:
            self.obj = lib_py_gel.Graph_copy(orig.obj)
    def __del__(self):
        lib_py_gel.Graph_delete(self.obj)
    def clear(self):
        """ Clear the graph. """
        lib_py_gel.Graph_clear(self.obj)
    def cleanup(self):
        """ Cleanup reorders the graph nodes such that there is no
        gap in the index range. """
        lib_py_gel.Graph_cleanup(self.obj)
    def nodes(self):
        """ Get all nodes as an iterable range """
        nodes = IntVector()
        lib_py_gel.Graph_nodes(self.obj, nodes.obj)
        return nodes
    def neighbors(self, n, mode='n'):
        """ Get the neighbors of node n. The final argument is either 'n' or 'e'. If it is 'n'
        the function returns all neighboring nodes, and if it is 'e' it returns incident edges."""
        nbors = IntVector()
        lib_py_gel.Graph_neighbors(self.obj, n, nbors.obj, ct.c_char(mode.encode('ascii')))
        return nbors
    def positions(self):
        """ Get the vertex positions by reference. You can assign to the
        positions. """
        pos = ct.POINTER(ct.c_double)()
        n = lib_py_gel.Graph_positions(self.obj, ct.byref(pos))
        return np.ctypeslib.as_array(pos,(n,3))
    def average_edge_length(self):
        """ Returns the average edge length. """
        ael = lib_py_gel.Graph_average_edge_length(self.obj)
        return ael
    def add_node(self, p):
        """ Adds node with position p to the graph and returns the
        index of the new node. """
        return lib_py_gel.Graph_add_node(self.obj, np.array(p))
    def remove_node(self, n):
        """ Removes the node n passed as argument. This does not change
        any indices of other nodes, but n is then invalid. """
        lib_py_gel.Graph_remove_node(self.obj, n)
    def node_in_use(self, n):
        """ Checks if n is in_use. This function returns false both
        if n has been removed and if n is an index outside the range of
        indices that are used. """
        return lib_py_gel.Graph_node_in_use(self.obj, n)
    def connect_nodes(self, n0, n1):
        """ Creates a new edge connecting nodes n0 and n1. The index of
        the new edge is returned. """
        return lib_py_gel.Graph_connect_nodes(self.obj, n0, n1)
    def disconnect_nodes(self, n0, n1):
        """ Disconect nodes n0 and n1"""
        lib_py_gel.Graph_disconnect_nodes(self.obj, n0, n1)
    def merge_nodes(self, n0, n1, avg_pos):
        """ Merge nodes n0 and n1. avg_pos indicates if you want the position to be the average. """
        lib_py_gel.Graph_merge_nodes(self.obj, n0, n1, avg_pos)


def from_mesh(m):
    """ Creates a graph from a mesh. The argument, m, is the input mesh,
    and the function returns a graph with the same vertices and edges
    as m."""
    g = Graph()
    lib_py_gel.graph_from_mesh(m.obj, g.obj)
    return g

def load(fn):
    """ Load a graph from a file. The argument, fn, is the filename which
    is in a special format similar to Wavefront obj. The loaded graph is
    returned by the function - or None if loading failed. """
    s = ct.c_char_p(fn.encode('utf-8'))
    g = Graph()
    if lib_py_gel.graph_load(g.obj, s):
        return g
    return None

def save(fn, g):
    """ Save graph to a file. The first argument, fn, is the file name,
    and g is the graph. This function returns True if saving happened and
    False otherwise. """
    s = ct.c_char_p(fn.encode('utf-8'))
    return lib_py_gel.graph_save(g.obj, s)

def to_mesh_cyl(g, fudge=0.0):
    """ Creates a Manifold mesh from the graph. The first argument, g, is the
    graph we want converted, and fudge is a constant that is used to increase the radius
    of every node. This is useful if the radii are 0. """
    m = hmesh.Manifold()
    lib_py_gel.graph_to_mesh_cyl(g.obj, m.obj, fudge)
    return m

def to_mesh_iso(g, fudge=0.0, res=256):
    """ Creates a Manifold mesh from the graph. The first argument, g, is the
    graph we want converted, and fudge is a constant that is used to increase the radius
    of every node. This is useful if the radii are 0. """
    m = hmesh.Manifold()
    lib_py_gel.graph_to_mesh_iso(g.obj, m.obj, fudge, res)
    return m


def smooth(g, iter=1, alpha=1.0):
    """ Simple Laplacian smoothing of a graph. The first argument is the Graph, g, iter
    is the number of iterations, and alpha is the weight. If the weight is high,
    each iteration causes a lot of smoothing, and a high number of iterations
    ensures that the effect of smoothing diffuses throughout the graph, i.e. that the
    effect is more global than local. """
    lib_py_gel.graph_smooth(g.obj, iter, alpha)

def edge_contract(g, dist_thresh):
    """ Simplifies a graph by contracting edges. The first argument, g, is the graph,
    and only edges shorter than dist_thresh are contracted. When an edge is contracted
    the merged vertices are moved to the average of their former positions. Thus,
    the ordering in which contractions are carried out matters. Hence, edges are
    contracted in the order of increasing length and edges are only considered if
    neither end point is the result of a contraction, but the process is then repeated
    until no more contractions are possible. Returns total number of contractions. """
    return lib_py_gel.graph_edge_contract(g.obj, dist_thresh)

def prune(g):
    """ Prune leaves of a graph. The graph, g, is passed as the argument. This function
        removes leaf nodes (valency 1) whose only neighbour has valency > 2. In practice
        such isolated leaves are frequently spurious if the graph is a skeleton. Does not
        return a value. """
    lib_py_gel.graph_prune(g.obj)
    
def saturate(g, hops=2, dist_frac=1.001, rad=1e300):
    """ Saturate the graph with edges. This is not a complete saturation. Edges are
    introduced between a vertex and other vertices that are reachable in hops steps, i.e.
    hops-order neighbors. dist_frac and rad are parameters used to govern the precise
    behaviour. Two nodes are only connected if their distance is less than rad and if
    their distance is less than dist_frac times the length of the path along existing
    edges in the graph. If dist_frac is at approximately 1 and rad is enormous, these
    two parameters make no difference. """
    lib_py_gel.graph_saturate(g.obj, hops, dist_frac, rad)
    
def LS_skeleton(g, sampling=True):
    """ Skeletonize a graph using the local separators approach. The first argument,
        g, is the graph, and, sampling indicates whether we try to use all vertices
        (False) as starting points for finding separators or just a sampling (True).
        The function returns a new graph which is the skeleton of the input graph. """
    skel = Graph()
    mapping = IntVector()
    lib_py_gel.graph_LS_skeleton(g.obj, skel.obj, mapping.obj, sampling)
    return skel
    
def LS_skeleton_and_map(g, sampling=True):
    """ Skeletonize a graph using the local separators approach. The first argument,
        g, is the graph, and, sampling indicates whether we try to use all vertices
        (False) as starting points for finding separators or just a sampling (True).
        The function returns a tuple containing a new graph which is the skeleton of
        the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    lib_py_gel.graph_LS_skeleton(g.obj, skel.obj, mapping.obj, sampling)
    return skel, mapping

def MSLS_skeleton(g, grow_thresh=64):
    """ Skeletonize a graph using the multi-scale local separators approach. The first argument,
        g, is the graph, and, sampling indicates whether we try to use all vertices
        (False) as starting points for finding separators or just a sampling (True).
        The function returns a new graph which is the skeleton of the input graph. """
    skel = Graph()
    mapping = IntVector()
    lib_py_gel.graph_MSLS_skeleton(g.obj, skel.obj, mapping.obj, grow_thresh)
    return skel
    
def MSLS_skeleton_and_map(g, grow_thresh=64):
    """ Skeletonize a graph using the multi-scale local separators approach. The first argument,
        g, is the graph, and, sampling indicates whether we try to use all vertices
        (False) as starting points for finding separators or just a sampling (True).
        The function returns a tuple containing a new graph which is the skeleton of
        the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    lib_py_gel.graph_MSLS_skeleton(g.obj, skel.obj, mapping.obj, grow_thresh)
    return skel, mapping


def front_skeleton_and_map(g, colors, intervals=100):
    """ Skeletonize a graph using the front separators approach. The first argument,
        g, is the graph, and, colors is an nD array where each column contains a sequence
        of floating point values - one for each node. We can have as many columns as needed
        for the front separator computation. We can think of this as a coloring
        of the nodes, hence the name. In practice, a coloring might just be the x-coordinate
        of the nodes or some other function that indicates something about the structure of the
        graph. The function returns a tuple containing a new graph which is the
        skeleton of the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    colors_flat = np.asarray(colors, dtype=ct.c_double, order='C')
    N_col = 1 if len(colors_flat.shape)==1 else colors_flat.shape[1]
    print("N_col:", N_col)
    pos = g.positions()
    lib_py_gel.graph_front_skeleton(g.obj, skel.obj, mapping.obj, N_col, colors_flat.ctypes.data_as(ct.POINTER(ct.c_double)), intervals)
    return skel, mapping

def front_skeleton(g, colors, intervals=100):
    """ Skeletonize a graph using the front separators approach. The first argument,
        g, is the graph, and, colors is a nD array where each column contains a sequence
        of floating point values - one for each node. We can have as many columns as needed
        for the front separator computation. We can think of this as a coloring
        of the nodes, hence the name. In practice, a coloring might just be the x-coordinate
        of the nodes or some other function that indicates something about the structure of the
        graph. The function returns a tuple containing a new graph which is the
        skeleton of the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    colors_flat = np.asarray(colors, dtype=ct.c_double, order='C')
    N_col = 1 if len(colors_flat.shape)==1 else colors_flat.shape[1]
    print("N_col:", N_col)
    pos = g.positions()
    lib_py_gel.graph_front_skeleton(g.obj, skel.obj, mapping.obj, N_col, colors_flat.ctypes.data_as(ct.POINTER(ct.c_double)), intervals)
    return skel

def combined_skeleton_and_map(g, colors, intervals=100):
    """ Skeletonize a graph using both the front separators approach and the multi scale local separators.
        The first argument, g, is the graph, and, colors is an nD array where each column contains a sequence
        of floating point values - one for each node. We can have as many columns as needed
        for the front separator computation. We can think of this as a coloring
        of the nodes, hence the name. In practice, a coloring might just be the x-coordinate
        of the nodes or some other function that indicates something about the structure of the
        graph. The function returns a tuple containing a new graph which is the
        skeleton of the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    colors_flat = np.asarray(colors, dtype=ct.c_double, order='C')
    N_col = 1 if len(colors_flat.shape)==1 else colors_flat.shape[1]
    print("N_col:", N_col)
    pos = g.positions()
    lib_py_gel.graph_combined_skeleton(g.obj, skel.obj, mapping.obj, N_col, colors_flat.ctypes.data_as(ct.POINTER(ct.c_double)), intervals)
    return skel, mapping

def combined_skeleton(g, colors, intervals=100):
    """ Skeletonize a graph using both the front separators approach and the multi scale local separators.
        The first argument, g, is the graph, and, colors is an nD array where each column contains a sequence
        of floating point values - one for each node. We can have as many columns as needed
        for the front separator computation. We can think of this as a coloring
        of the nodes, hence the name. In practice, a coloring might just be the x-coordinate
        of the nodes or some other function that indicates something about the structure of the
        graph. The function returns a new graph which is the
        skeleton of the input graph and a map from the graph nodes to the skeletal nodes. """
    skel = Graph()
    mapping = IntVector()
    colors_flat = np.asarray(colors, dtype=ct.c_double, order='C')
    N_col = 1 if len(colors_flat.shape)==1 else colors_flat.shape[1]
    print("N_col:", N_col)
    pos = g.positions()
    lib_py_gel.graph_combined_skeleton(g.obj, skel.obj, mapping.obj, N_col, colors_flat.ctypes.data_as(ct.POINTER(ct.c_double)), intervals)
    return skel
