"""Radius-limited Delaunay or Gabriel edges.

    python delaunay_edges.py [--gabriel] [mesh] [radius_factor]

With no mesh, 200 points on the unit sphere are used. The radius is 10,
unless a factor is given after the mesh path. With a mesh, the radius is
`radius_factor` times the average edge length, and the factor defaults to 4.

`--gabriel` keeps the edges whose diametral ball contains no other point.
"""

import sys

import numpy as np
from pygel3d import gl_display as gl
from pygel3d import graph, hmesh


def sphere_points(n=200, seed=1):
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    return points


def edge_count(g):
    count = 0
    for node in g.nodes():
        for neighbor in g.neighbors(node):
            if node < neighbor:
                count += 1
    return count


def main():
    gabriel = "--gabriel" in sys.argv
    args = [arg for arg in sys.argv[1:] if arg != "--gabriel"]
    if args:
        mesh = hmesh.load(args[0])
        if mesh is None:
            raise SystemExit(f"could not load {args[0]}")
        points = np.array(mesh.positions(), copy=True)
        factor = float(args[1]) if len(args) > 1 else 4.0
        radius = factor * hmesh.average_edge_length(mesh)
    else:
        points = sphere_points()
        radius = 10.0

    if gabriel:
        g = graph.gabriel_edges(points, radius)
        kind = "Gabriel"
    else:
        g = graph.delaunay_edges(points, radius)
        kind = "Delaunay"
    print(f"{points.shape[0]} points, radius {radius:.6g}, {edge_count(g)} {kind} edges")
    gl.Viewer().display(g)


if __name__ == "__main__":
    main()
