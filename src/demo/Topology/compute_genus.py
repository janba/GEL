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

components = hmesh.analyze_topology(m)

for i,comp in enumerate(components):
    print(f"Component {i+1}:")
    print(f"- Vertices: {comp['V']}")
    print(f"- Edges: {comp['E']}")
    print(f"- Faces: {comp['F']}")
    print(f"- Number of boundary curves: {comp['b']}")
    print(f"- Genus: {comp['g']}")