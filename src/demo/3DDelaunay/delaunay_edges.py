"""Radius-limited Delaunay or Gabriel edges.

    python delaunay_edges.py [--gabriel] [--mst] [--no-tunnel-weight] [--radius FACTOR] [mesh] [radius_factor]

With no mesh, 200 points on the unit sphere are used and the radius is 10.
`--radius` is the factor on length, default 4. A mesh uses that factor times
the average edge length. A file whose mesh has no vertices is read as a point
cloud: `.obj` uses each line whose first word is `v`, and any other file uses
the first three numbers on each line. The radius is then the factor times the
median distance to the nearest other point. A factor after the mesh path is
the same number; `--radius` wins when both are given.

`--gabriel` keeps the edges whose diametral ball contains no other point.
`--mst` draws a minimum spanning forest of that graph. The forest and the
tunnel lengths are computed in C++. Each edge length is
`length * (1 + tunnel_weight)`. The weight is the RsR tunnel alignment on
the unit edge, `|ê·n| + |ê·m|`. RsR divides that by the neighborhood radius
when the term is added to an angle; the length is already outside here, so
an edge along both normals is three times its Euclidean length and a tangent
edge is unchanged. `--no-tunnel-weight` leaves the length Euclidean. Normals
are the adaptive covariance estimate with the covariance sign left in place.
Each executed step prints its time. Point loading is not timed.
"""

import re
import sys
import time
from os.path import splitext

import numpy as np
from pygel3d import gl_display as gl
from pygel3d import graph, hmesh

# A number, and the text around it is the separator. Extra numbers are ignored.
_NUMBER = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?")


def sphere_points(n=200, seed=1):
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    return points


def parse_args(argv):
    flags = {"--gabriel", "--mst", "--no-tunnel-weight"}
    selected = set()
    positional = []
    radius_factor = None
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg == "--radius":
            if index + 1 >= len(argv):
                raise SystemExit("--radius needs a factor")
            radius_factor = float(argv[index + 1])
            index += 2
            continue
        if arg in flags:
            selected.add(arg)
        else:
            positional.append(arg)
        index += 1
    return selected, positional, radius_factor


def read_point_cloud(path):
    """Points from a file that did not load as a mesh with vertices.

    `.obj` keeps lines whose first word is `v`, then the next three
    whitespace-separated coordinates. Any other extension takes the first
    three numbers on each line. A separator is whatever is not itself a
    number. Lines with fewer than three coordinates are skipped.
    """
    obj = splitext(path)[1].lower() == ".obj"
    points = []
    with open(path) as source:
        for line in source:
            if obj:
                parts = line.split()
                if len(parts) < 4 or parts[0] != "v":
                    continue
                coords = parts[1:4]
            else:
                coords = _NUMBER.findall(line)[:3]
                if len(coords) < 3:
                    continue
            points.append(tuple(float(value) for value in coords))
    if not points:
        return np.zeros((0, 3))
    return np.asarray(points, dtype=np.float64)


def neighbor_spacing(points):
    """Median distance from each point to the nearest other point."""
    if len(points) < 2:
        return 1.0
    from scipy.spatial import KDTree

    distances, _ = KDTree(points).query(points, k=2)
    nearest = np.asarray(distances[:, 1], dtype=np.float64)
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if nearest.size == 0:
        return 1.0
    return float(np.median(nearest))


def load_points(positional, radius_factor=None):
    if not positional:
        return sphere_points(), 10.0
    path = positional[0]
    if radius_factor is None:
        factor = float(positional[1]) if len(positional) > 1 else 4.0
    else:
        factor = radius_factor
    mesh = hmesh.load(path)
    if mesh is not None and len(mesh.vertices()) > 0:
        points = np.array(mesh.positions(), copy=True)
        return points, factor * hmesh.average_edge_length(mesh)
    points = read_point_cloud(path)
    if len(points) == 0:
        raise SystemExit(f"could not load {path}")
    return points, factor * neighbor_spacing(points)


def format_duration(seconds):
    if seconds < 0.001:
        return f"{seconds * 1e6:.0f} us"
    if seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    return f"{seconds:.2f} s"


def edges_to_show(points, radius, gabriel, show_mst, tunnel):
    """Build the graph to draw and print each step that runs."""
    started = time.perf_counter()
    if gabriel:
        g = graph.gabriel_edges(points, radius)
        kind = "Gabriel"
    else:
        g = graph.delaunay_edges(points, radius)
        kind = "Delaunay"
    print(
        f"{points.shape[0]} points, radius {radius:.6g}, {g.no_edges()} {kind} edges "
        f"({format_duration(time.perf_counter() - started)})"
    )

    shown = g
    if show_mst:
        shown, weights, normal_seconds, forest_seconds = graph.minimum_spanning_forest(g, tunnel)
        if normal_seconds >= 0:
            print(f"normal estimation {format_duration(normal_seconds)}")
        if tunnel and g.no_edges():
            print(
                "tunnel weight "
                f"min {weights[0]:.3g}, median {weights[1]:.3g}, max {weights[2]:.3g}"
            )
        forest_edges = shown.no_edges()
        components = points.shape[0] - forest_edges
        extra = f", {components} components" if components > 1 else ""
        description = "tunnel-weighted length" if tunnel else "Euclidean length"
        print(f"MST {forest_edges} edges, {description}{extra} ({format_duration(forest_seconds)})")
    return shown


def main():
    selected, positional, radius_factor = parse_args(sys.argv[1:])
    points, radius = load_points(positional, radius_factor)
    shown = edges_to_show(
        points,
        radius,
        "--gabriel" in selected,
        "--mst" in selected,
        "--no-tunnel-weight" not in selected,
    )
    gl.Viewer().display(shown)


if __name__ == "__main__":
    main()
