""" This is a module with a function, display, that provides functionality for displaying a
    Manifold or a Graph as an interactive 3D model in a Jupyter Notebook. """
from pygel3d import hmesh, graph
from numpy import array
import plotly.graph_objs as go
import plotly.offline as py

EXPORT_MODE = False

def set_export_mode(_exp_mode=True):
    """ Calling this function will set export mode to true. It is necessary
    to do so if we wish to export a notebook containing interactive
    plotly graphics (made with display below) to HTML. In other words, this function
    should not necessarily be called in normal usage but only when we export to HTML. It is
    then called once in the beginning of the notebook. However, as a bit of a twist on
    this story, it appears that if we don't call this function, any call to display must
    be the last thing that happens in a cell. So, maybe it is best to always call
    set_export_mode in the beginning of a notebook.
    """
    global EXPORT_MODE
    EXPORT_MODE=_exp_mode
    if EXPORT_MODE:
        py.init_notebook_mode(connected=False)

def display(m,wireframe=True,smooth=True,data=None):
    """ The display function shows an interactive presentation of the Manifold, m, inside
        a Jupyter Notebook. wireframe=True means that a wireframe view of the mesh is
        superimposed on the 3D model. If smooth=True, the mesh is rendered with vertex
        normals. Otherwise, the mesh is rendered with face normals. If data=None, the
        mesh is shown in a light grey color. If data contains an array of scalar values
        per vertex, these are mapped to colors used to color the mesh. Finally, note that
        m can also be a Graph. In that case the display function just draws the edges as
        black lines. """
    mesh_data = []
    if isinstance(m,hmesh.Manifold):
        xyz = array([ p for p in m.positions()])
        m_tri = hmesh.Manifold(m)
        hmesh.triangulate(m_tri)
        ijk = array([[ idx for idx in m_tri.circulate_face(f,'v')] for f in m_tri.faces()])
        mesh = go.Mesh3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],
                i=ijk[:,0],j=ijk[:,1],k=ijk[:,2],color='#dddddd',flatshading=not smooth)
        if data is not None:
            mesh['intensity'] = data
        mesh_data += [mesh]
        if wireframe:
            pos = m.positions()
            xyze = []
            for h in m.halfedges():
                if h < m.opposite_halfedge(h):
                    p0 = pos[m.incident_vertex(m.opposite_halfedge(h))]
                    p1 = pos[m.incident_vertex(h)]
                    xyze.append(array(p0))
                    xyze.append(array(p1))
                    xyze.append(array([None, None, None]))
            xyze = array(xyze)
            trace1=go.Scatter3d(x=xyze[:,0],y=xyze[:,1],z=xyze[:,2],
                       mode='lines',
                       line=dict(color='rgb(125,0,0)', width=1),
                       hoverinfo='none')
            mesh_data += [trace1]
    elif isinstance(m,graph.Graph):
        pos = m.positions()
        xyze = []
        for v in m.nodes():
            for w in m.neighbors(v):
                if v < w:
                    p0 = pos[v]
                    p1 = pos[w]
                    xyze.append(array(p0))
                    xyze.append(array(p1))
                    xyze.append(array([None, None, None]))
        xyze = array(xyze)
        trace1=go.Scatter3d(x=xyze[:,0],y=xyze[:,1],z=xyze[:,2],
                   mode='lines',
                   line=dict(color='rgb(0,0,0)', width=1),
                   hoverinfo='none')
        mesh_data += [trace1]

        
    lyt = go.Layout(width=850,height=800)
    lyt.scene.aspectmode="data"
    if EXPORT_MODE:
        py.iplot(dict(data=mesh_data,layout=lyt))
    else:
        return go.FigureWidget(mesh_data,lyt)

