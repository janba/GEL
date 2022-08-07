""" This modules provides an OpenGL based viewer for graphs and meshes """
from pygel3d import hmesh, graph, lib_py_gel
import ctypes as ct
import numpy as np
from os import getcwd, chdir

try:
    lib_py_gel.GLManifoldViewer_new.restype = ct.c_void_p
    lib_py_gel.GLManifoldViewer_delete.argtypes = (ct.c_void_p,)
    lib_py_gel.GLManifoldViewer_display.argtypes = (ct.c_void_p,ct.c_void_p,ct.c_void_p,ct.c_char,ct.c_bool, ct.POINTER(ct.c_float*3), ct.POINTER(ct.c_double),ct.c_bool,ct.c_bool)
    lib_py_gel.GLManifoldViewer_get_annotation_points.restype = ct.c_size_t
    lib_py_gel.GLManifoldViewer_get_annotation_points.argtypes = (ct.c_void_p, ct.POINTER(ct.POINTER(ct.c_double)))
    lib_py_gel.GLManifoldViewer_set_annotation_points.argtypes = (ct.c_void_p, ct.c_int, ct.POINTER(ct.c_double))
    lib_py_gel.GLManifoldViewer_event_loop.argtypes = (ct.c_bool,)
    class Viewer:
        """ An OpenGL Viewer for Manifolds and Graphs. Having created an instance of this
        class, call display to show a mesh or a graph. The display function is flexible,
        allowing several types of interactive visualization. Each instance of this
        class corresponds to a single window, but you can have several
        GLManifoldViewer and hence also several windows showing different
        visualizations. """
        def __init__(self):
            current_directory = getcwd()
            self.obj = lib_py_gel.GLManifoldViewer_new()
            chdir(current_directory) # Necessary because init_glfw changes cwd
        def __del__(self):
            lib_py_gel.GLManifoldViewer_delete(self.obj)
        def display(self, m, g=None, mode='w', smooth=True, bg_col=[0.3,0.3,0.3], data=None, reset_view=False, once=False):
            """ Display a mesh

            Args:
            ---
            - m : the Manifold mesh or Graph we want to show.
            - g : the Graph we want to show. If you only want to show a graph, you
                can simply pass the graph as m, so the g argument is relevant only if
                you need to show both a Manifold _and_ a Graph.
            - mode : a single character that determines how the mesh is visualized:
                'w' - wireframe,
                'i' - isophote,
                'g' - glazed (try it and see),
                's' - scalar field,
                'l' - line field,
                'n' - normal.
                'x' - xray or ghost rendering. Useful to show Manifold on top of Graph
            - smooth : if True we use vertex normals. Otherwise, face normals.
            - bg_col : background color.
            - data : per vertex data for visualization. scalar or vector field.
            - reset_view : if False view is as left in the previous display call. If
                True, the view is reset to the default.
            - once : if True we immediately exit the event loop and return. However,
                the window stays and if the event loop is called from this or any
                other viewer, the window will still be responsive.

            Interactive controls:
            ---
            When a viewer window is displayed on the screen, you can naviagate with
            the mouse: Left mouse button rotates, right mouse button is used for
            zooming and (if shift is pressed) for panning. If you hold control, any
            mouse button will pick a point on the 3D model. Up to 19 of these points
            have unique colors.  If you pick an already placed annotation point it
            will be removed and can now be placed elsewhere. Hit space bar to clear
            the annotation points. Hitting ESC exits the event loop causing control
            to return to the script.
            """
            data_ct = np.array(data,dtype=ct.c_double).ctypes
            data_a = data_ct.data_as(ct.POINTER(ct.c_double))
            bg_col_ct = np.array(bg_col,dtype=ct.c_float).ctypes
            bg_col_a = bg_col_ct.data_as(ct.POINTER(ct.c_float*3))
            if isinstance(m,graph.Graph):
                g = m
                m = None
            if isinstance(m,hmesh.Manifold) and isinstance(g,graph.Graph):
                lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, g.obj, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
            elif isinstance(m,hmesh.Manifold):
                lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, 0, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
            elif isinstance(g,graph.Graph):
                lib_py_gel.GLManifoldViewer_display(self.obj, 0, g.obj, ct.c_char(mode.encode('ascii')),smooth,bg_col_a,data_a,reset_view,once)
                
        def annotation_points(self):
            """ Retrieve a vector of annotation points. This vector is not a copy,
            so any changes made to the points will be reflected in the viewer. """
            pos = ct.POINTER(ct.c_double)()
            n = lib_py_gel.GLManifoldViewer_get_annotation_points(self.obj, ct.byref(pos))
            if n == 0:
                return None
            return np.ctypeslib.as_array(pos,(n,3))
        def set_annotation_points(self, pts):
            n = int(np.size(pts)/3)
            pts_ct = np.array(pts,dtype=ct.c_double).ctypes
            pts_a = pts_ct.data_as(ct.POINTER(ct.c_double))
            lib_py_gel.GLManifoldViewer_set_annotation_points(self.obj, n, pts_a)
        @staticmethod
        def event_loop():
            """ Explicit call to the event loop. This function enters the event loop.
            Call it if you want to turn on interactivity in the currently displayed
            window."""
            lib_py_gel.GLManifoldViewer_event_loop(False)
except AttributeError:
    pass
