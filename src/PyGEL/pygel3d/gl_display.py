""" This modules provides an OpenGL based viewer for graphs and meshes """
from pygel3d import lib_py_gel
from pygel3d.hmesh import Manifold
from pygel3d.graph import Graph
import numpy as np
from os import getcwd, chdir

try:
    # GL viewer functionality with pybind11
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

        def clone_controller(self, v):
            """ Clone another viewers navigation controller. This function
            can be used to synchronize the view between several viewers by calling
            this function on all but the first viewer after each call to display.
            It works by giving several viewer objects the same trackball controller. """
            lib_py_gel.GLManifoldViewer_clone_controller(self.obj, v.obj)

        def display(self, m=Manifold(), g=Graph(), mode='w', smooth=True, data=np.zeros(1), bg_col=[0,0,0], reset_view=True, once=True):
            """ The swiss army knife of display methods. This method displays
            both Manifolds and Graphs. Each argument is optional. This function
            is very flexible and hence slightly complex:

            m is a Manifold. If provided it will be displayed.

            g is a Graph. If provided it will be displayed.

            mode specifies how mesh or graph is displayed. The following modes are
            available
            'w' - wireframe mode.
            's' - shaded mode using silhouette like, lit rendering.
            'n' - show the mesh as points with normal vectors
            'c' - color the mesh according to discrete Gaussian curvature.

            smooth governs whether the mesh is smoothed before display. In some
            cases smoothed display gives a better impression of the shape, but
            it is not representative of the actual geometry.

            data is a scalar per vertex data vector. If provided, vertices are
            color coded accoding to the data values.

            bg_col specifies the background color of the scene. default=[0,0,0]

            reset_view specifies whether to reset the view when the object is
            first displayed. default=True. In general it is safe to leave this
            setting as is, and if you call the display function on the same Viewer
            instance with the same (possibly modified) model, it is usually possible
            to inspect the modificaitons without the view changes. Sometimes, however
            it might be useful to be able to prevent the view from being reset.

            once controls whether the function returns immediately or after the
            window has been closed by the user. If True, the function returns immediately.
            If False, the viewer enters an event loop and will return when the user closes
            the window. default=True. """
            lib_py_gel.GLManifoldViewer_display(self.obj, m.obj)
            # if len(data) > 1:
            #     data_np = np.array(data, dtype=np.float64)
            #     bg_col_np = np.array(bg_col, dtype=np.float32) if len(bg_col) == 3 else np.array([0,0,0], dtype=np.float32)
            #     lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, g.obj, mode.encode('ascii'), smooth, bg_col_np, data_np, reset_view, once)
            # else:
            #     bg_col_np = np.array(bg_col, dtype=np.float32) if len(bg_col) == 3 else np.array([0,0,0], dtype=np.float32)
            #     lib_py_gel.GLManifoldViewer_display(self.obj, m.obj, g.obj, mode.encode('ascii'), smooth, bg_col_np, None, reset_view, once)

        def get_annotation_points(self):
            """ Get the annotation points. It is possible to add points to the
            display by shift-clicking in the viewer. These points are returned
            by this function as a numpy array. """
            points = lib_py_gel.GLManifoldViewer_get_annotation_points(self.obj)
            return np.array(points)

        def set_annotation_points(self, pts):
            """ Set the annotation points. These can be shown by pressing 'p' in
            the viewer. """
            pts_np = np.array(pts, dtype=np.float64)
            lib_py_gel.GLManifoldViewer_set_annotation_points(self.obj, pts_np)

    # Make the event loop function available at the module level
    def event_loop(once=False):
        """ Run the GLFW event loop. If once is True, the event loop polls events
        and returns immediately. Otherwise, the function waits for events. """
        lib_py_gel.GLManifoldViewer_event_loop(once)

except ImportError:
    # Fallback when OpenGL viewer is not available
    print("OpenGL viewer not available. Please install dependencies.")
    
    class Viewer:
        def __init__(self):
            print("OpenGL viewer not available")
        def display(self, *args, **kwargs):
            print("OpenGL display not available")
        def get_annotation_points(self):
            return np.array([])
        def set_annotation_points(self, pts):
            pass
    
    def event_loop(once=False):
        print("OpenGL event loop not available")
