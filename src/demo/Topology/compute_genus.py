#!/opt/local/bin/python
from sys import argv
from numpy import zeros
from pygel3d import hmesh

def connected_components(m: hmesh.Manifold):
    """
    Segments the mesh into connected components based on face traversal.
    For each component, creates a new separate mesh and returns a list of meshes.
    """
    visited = zeros(m.no_allocated_faces(), dtype=bool)
    components = []
    pos = m.positions()
    for face in m.faces():
        if not visited[face]:
            # Start a new component
            stack = [face]
            # Create a new mesh for the component
            component_mesh = hmesh.Manifold()
            while stack:
                current_face = stack.pop()
                if not visited[current_face]:
                    visited[current_face] = True
                    component_mesh.add_face([ pos[v] for v in m.circulate_face(current_face) ])
                    # Add neighboring faces to the stack
                    for neighbor in m.circulate_face(current_face, mode='f'):
                        if m.face_in_use(neighbor) and not visited[neighbor]:
                            stack.append(neighbor)
            hmesh.stitch(component_mesh)
            components.append(component_mesh)

    return components

def boundary_curves(m: hmesh.Manifold):
    """
    Extracts the boundary curves of the mesh.
    Returns a list of boundary curves, each represented as a list of vertices.
    """
    boundary_curves = []
    visited = zeros(m.no_allocated_halfedges(), dtype=bool)
    for h0 in m.halfedges():
        if visited[h0]:
            continue
        if not m.face_in_use(m.incident_face(h0)):
            bc = [] # Boundary curve
            h = h0
            while not visited[h]:
                visited[h] = True
                bc.append(h)
                h = m.next_halfedge(h)
            boundary_curves.append(bc)
    return boundary_curves


if __name__ == "__main__":
    # Load the mesh file
    mesh_file = argv[1]
    m = hmesh.load(mesh_file)

    components = hmesh.connected_components(m)
    
    print(f"Number of connected components: {len(components)}")


    for i,comp in enumerate(components):
        print(f"Component {i+1}:")
        b = len(hmesh.count_boundary_curves(comp))
        V = len(comp.vertices())  # Number of vertices
        E = len(comp.halfedges())//2  # Number of edges
        F = len(comp.faces())  # Number of faces
        print(f"Vertices: {V}")
        print(f"Edges: {E}")
        print(f"Faces: {F}")
        print(f"Number of boundary curves: {b}")
        # V - E + F = 2 - 2g - b
        # g = -(V - E + F - 2 + b)/2
        g = -(V - E + F - 2 + b)//2
        print(f"Genus: {g}")