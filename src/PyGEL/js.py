""" PyGEL.js is a module with a single function, display, that provides functionality for displaying a mesh
    Manifold as an interactive 3D model in a Jupyter Notebook """
from PyGEL import gel
from numpy import array
import plotly.offline as py
import plotly.graph_objs as go
py.init_notebook_mode(connected=True)

def display(m,wireframe=True,smooth=True,data=None):
    """ The display function shows an interactive presentation of the Manifold, m, inside
        a Jupyter Notebook. wireframe=True means that a wireframe view of the mesh is
        superimposed on the 3D model. If smooth=True, the mesh is rendered with vertex
        normals. Otherwise, the mesh is rendered with face normals. If data=None, the
        mesh is shown in a light grey color. If data contains an array of scalar values
        per vertex, these are mapped to colors used to color the mesh."""
    xyz = array([ p for p in m.positions()])
    m_tri = gel.Manifold(m)
    gel.triangulate(m_tri)
    ijk = array([[ idx for idx in m_tri.circulate_face(f,'v')] for f in m_tri.faces()])
    mesh = go.Mesh3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],
            i=ijk[:,0],j=ijk[:,1],k=ijk[:,2],color='#dddddd',flatshading=not smooth)
    if data is not None:
        mesh['intensity'] = data
    mesh_data = [mesh]
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
                   line=go.Line(color='rgb(125,0,0)', width=1),
                   hoverinfo='none')
        mesh_data += [trace1]
    lyt = go.Layout(scene=go.Scene(aspectmode='data'))
    fig = go.Figure(data=mesh_data,layout=lyt)
    py.iplot(fig)
