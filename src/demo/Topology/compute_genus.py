#!/opt/local/bin/python
# This very simple script loads a mesh and splits it into its connected components.
# It then computes the number of vertices, edges, faces, boundary curves, and genus
# for each component. 
#
# Usage: python compute_genus.py <mesh_file>

from sys import argv
from numpy import zeros
from pygel3d import hmesh

# Load the mesh file
m = hmesh.load(argv[1])
if m is None:
    print("Failed to load mesh.")
    exit(1)

components = hmesh.connected_components(m)

print(f"Number of connected components: {len(components)}")
for i,comp in enumerate(components):
    print(f"Component {i+1}:")
    b = hmesh.count_boundary_curves(comp)
    V = len(comp.vertices())  # Number of vertices
    E = len(comp.halfedges())//2  # Number of edges
    F = len(comp.faces())  # Number of faces
    g = -(V - E + F - 2 + b)//2 # Genus calculation
    print(f"Vertices: {V}")
    print(f"Edges: {E}")
    print(f"Faces: {F}")
    print(f"Number of boundary curves: {b}")
    print(f"Genus: {g}")