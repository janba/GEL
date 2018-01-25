#from pythreejs import *
#from IPython.display import display
#from ipywidgets import HTML, Text
#from traitlets import link, dlink
import PyGEL as pgl
from numpy import array
import plotly.offline as py
import plotly.graph_objs as go
py.init_notebook_mode(connected=False)


def display_js(m,wireframe=True,smooth=True,data=None):
#    pmin,pmax = pgl.bbox(m)
#    aspect = list(pmax-pmin)
#    aspect /= aspect[0]
    xyz = array([ p for p in m.positions()])
    ijk = array([[ idx for idx in m.circulate_face(f,'v')] for f in m.faces()])
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
