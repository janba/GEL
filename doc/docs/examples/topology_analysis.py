# Usage: python topology_analysis.py <input.obj>
import pygel3d.hmesh as hmesh
from sys import argv

def analyze_mesh(filename):
    m = hmesh.load(filename)
    
    # Topology
    components = hmesh.analyze_topology(m)
    for i,comp in enumerate(components):
        print(f"Component {i+1}:")
        print(f"- Vertices: {comp['V']}")
        print(f"- Edges: {comp['E']}")
        print(f"- Faces: {comp['F']}")
        print(f"- Number of boundary curves: {comp['b']}")
        print(f"- Genus: {comp['g']}")

    # Geometry
    bbox_min, bbox_max = hmesh.bbox(m)
    print(f"  Bounding box: {bbox_min} to {bbox_max}")
    
    sphere_center, sphere_radius = hmesh.bsphere(m)
    print(f"  Bounding sphere: center {sphere_center}, radius {sphere_radius}")
    
    print(f"  Total area: {hmesh.area(m):.4f}")
    if hmesh.closed(m):
        print(f"  Volume: {hmesh.volume(m):.4f}")
    
    # Edge statistics
    avg_len = hmesh.average_edge_length(m)
    med_len = hmesh.median_edge_length(m)
    print(f"  Average edge length: {avg_len:.4f}")
    print(f"  Median edge length: {med_len:.4f}")

# Analyze a mesh
analyze_mesh(argv[1])
