#from pythreejs import *
#from IPython.display import display
#from ipywidgets import HTML, Text
#from traitlets import link, dlink
import PyGEL as pgl
from numpy import array
import plotly.plotly as py
import plotly
import plotly.graph_objs as go

import plotly.offline as py
import plotly.graph_objs as go
py.init_notebook_mode(connected=False)


def display_js(m,wireframe=False,flatshading=True):
#    pmin,pmax = pgl.bbox(m)
#    aspect = list(pmax-pmin)
#    aspect /= aspect[0]
    xyz = array([ p for p in m.positions()])
    ijk = array([[ idx for idx in m.circulate_face(f,'v')] for f in m.faces()])
    mesh = go.Mesh3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],
            i=ijk[:,0],j=ijk[:,1],k=ijk[:,2], intensity=xyz[:,2],flatshading=flatshading)
    data = [mesh]
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
        data += [trace1]
    lyt = go.Layout(scene=go.Scene(aspectmode='data'))
    fig = go.Figure(data=data,layout=lyt)
    py.iplot(fig)
# def display_js(m,wireframe=True,vertex_normals=False,color=0x557799):
#     verts = [ coord for p in m.positions() for coord in list(p) ]
#     faces = [ v for f in m.faces() for v in m.circulate_face(f,'v')]
#     plot = k3d.plot()
#     mesh = k3d.mesh(verts, faces, color)
#     plot += mesh
#     plot.display()

# def display_js(m,wireframe=True,vertex_normals=False,color='#eebb99'):
#     c,r = pgl.bsphere(m)
#     r = 2.0 * r.value
#     pos = [-c[0],-c[1],-c[2]]
#     verts = [ coord for p in m.positions() for coord in list(p) ]
#     faces = [ v for f in m.faces() for v in m.circulate_face(f,'v')]
#     geo = FaceGeometry(vertices=verts,face3=faces,faceNormals=True)
#     m_solid = Mesh(geometry=geo,material=PhongMaterial(color=color),position=pos)
#     scene_objects = [AmbientLight(color='#777777'), m_solid]
#     if wireframe:
#         m_wire = Mesh(geometry=geo, material=BasicMaterial(color='#000033',wireframe=True), position=pos)
#         scene_objects +=  [m_wire]
#     scene = Scene(children=scene_objects)
#     c = PerspectiveCamera(position=[0, 0, r], up=[0, 0, 1],
#                       children=[DirectionalLight(color='white',
#                                                  position=[0, 0, r],
#                                                  intensity=0.5)])
#     renderer = Renderer(camera=c,
#                     scene=scene,
#                     controls=[OrbitControls(controlling=c)],
#                     antialias=True)
#     display(renderer)
