from pythreejs import *
from IPython.display import display
from ipywidgets import HTML, Text
from traitlets import link, dlink
import PyGEL as pgl
from numpy import array

def display_js(m,wireframe=True,vertex_normals=False,color='#eebb99'):
    c,r = pgl.bsphere(m)
    r = 2.0 * r.value
    pos = [-c[0],-c[1],-c[2]]
    verts = [ coord for p in m.positions() for coord in list(p) ]
    faces = [ v for f in m.faces() for v in m.circulate_face(f,'v')]
    geo = FaceGeometry(vertices=verts,face3=faces,faceNormals=True)
    m_solid = Mesh(geometry=geo,material=PhongMaterial(color=color),position=pos)
    scene_objects = [AmbientLight(color='#777777'), m_solid]
    if wireframe:
        m_wire = Mesh(geometry=geo, material=BasicMaterial(color='#000033',wireframe=True), position=pos)
        scene_objects +=  [m_wire]
    scene = Scene(children=scene_objects)
    c = PerspectiveCamera(position=[0, 0, r], up=[0, 0, 1],
                      children=[DirectionalLight(color='white',
                                                 position=[0, 0, r],
                                                 intensity=0.5)])
    renderer = Renderer(camera=c,
                    scene=scene,
                    controls=[OrbitControls(controlling=c)],
                    antialias=True)
    display(renderer)
