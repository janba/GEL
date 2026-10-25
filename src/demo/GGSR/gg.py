"""Visualize the radius-limited Gabriel graph.

    python gg.py [--radius FACTOR] [points.obj]

With no file, 400 points on the unit sphere are shown. The radius is
FACTOR times the median nearest-neighbor distance. The default factor
is 4. An .obj point cloud uses each line whose first word is v.
"""

import sys
from os.path import splitext

import numpy as np
import polyscope as ps
from pygel3d.graph import gabriel_edges


def sphere_points(n=400, seed=1):
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    return points


def read_obj_vertices(path):
    points = []
    with open(path) as source:
        for line in source:
            parts = line.split()
            if len(parts) < 4 or parts[0] != "v":
                continue
            points.append(tuple(float(value) for value in parts[1:4]))
    if not points:
        return np.zeros((0, 3))
    return np.asarray(points, dtype=np.float64)


def parse_args(argv):
    radius_factor = 4.0
    positional = []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg == "--radius":
            if index + 1 >= len(argv):
                raise SystemExit("--radius needs a factor")
            radius_factor = float(argv[index + 1])
            index += 2
            continue
        positional.append(arg)
        index += 1
    return positional, radius_factor


def median_neighbor_spacing(points):
    from scipy.spatial import KDTree

    distances, _ = KDTree(points).query(points, k=2)
    nearest = np.asarray(distances[:, 1], dtype=np.float64)
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if nearest.size == 0:
        return 1.0
    return float(np.median(nearest))


def gabriel_segments(points, radius):
    """Undirected Gabriel edges as an (E, 2) array of point indices."""
    graph = gabriel_edges(points, radius)
    segments = []
    for index in range(len(points)):
        for neighbor in graph.neighbors(index):
            other = int(neighbor)
            if index < other:
                segments.append((index, other))
    return graph, np.asarray(segments, dtype=np.int32).reshape(-1, 2)


def main():
    positional, radius_factor = parse_args(sys.argv[1:])
    if not positional:
        points = sphere_points()
        source = "unit sphere"
    else:
        path = positional[0]
        if splitext(path)[1].lower() != ".obj":
            raise SystemExit(f"expected an .obj point cloud, got {path}")
        points = read_obj_vertices(path)
        if len(points) == 0:
            raise SystemExit(f"could not load {path}")
        source = path

    radius = radius_factor * median_neighbor_spacing(points)
    graph, segments = gabriel_segments(points, radius)
    print(
        f"{source}: {len(points)} points, radius factor {radius_factor}, "
        f"radius {radius:.6g}, {graph.no_edges()} Gabriel edges"
    )

    ps.init()
    ps.set_ground_plane_mode("shadow_only")
    ps.register_point_cloud("points", points, radius=0.004)
    if len(segments):
        ps.register_curve_network("Gabriel graph", points, segments, radius=0.0015)
    ps.show()


if __name__ == "__main__":
    main()
