"""Visualize the radius-limited Gabriel graph.

    python gg.py [--radius FACTOR] [--remove-tunnel [ALIGNMENT]] [points.obj]

With no file, 400 points on the unit sphere are shown. The radius is
FACTOR times the median nearest-neighbor distance. The default factor
is 4. An .obj point cloud uses each line whose first word is v.

`--remove-tunnel` colors an edge red when it leaves the tangent plane at
both ends. `|n·ê|` and `|m·ê|` must both be at least ALIGNMENT. The
default is 0.5. Normals are the unoriented Euclidean covariance
estimates used by GGSR. The other Gabriel edges keep their color.
"""

import sys
from os.path import splitext

import numpy as np
import polyscope as ps
from pygel3d.graph import gabriel_edges
from pygel3d.hmesh import rsr_normals


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


def float_token(text):
    try:
        value = float(text)
    except ValueError:
        return None
    if not np.isfinite(value):
        return None
    return value


def parse_args(argv):
    radius_factor = 4.0
    tunnel_alignment = None
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
        if arg == "--remove-tunnel":
            tunnel_alignment = 0.5
            if index + 1 < len(argv):
                value = float_token(argv[index + 1])
                if value is not None:
                    tunnel_alignment = value
                    index += 1
            if not 0.0 < tunnel_alignment <= 1.0:
                raise SystemExit("--remove-tunnel needs an alignment in (0, 1]")
            index += 1
            continue
        positional.append(arg)
        index += 1
    return positional, radius_factor, tunnel_alignment


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


def tunnel_edges(points, segments, alignment):
    """True where the edge leaves the tangent plane at both ends.

    `|n·ê|` and `|m·ê|` must both be at least `alignment`. Normals are the
    unoriented Euclidean covariance estimates used by GGSR. A zero normal
    does not mark its edges.
    """
    mask = np.zeros(len(segments), dtype=bool)
    if len(segments) == 0:
        return mask
    normals = np.asarray(
        rsr_normals(points, use_Euclid_dist=True, orientation="none"),
        dtype=np.float64,
    )
    ends0 = segments[:, 0]
    ends1 = segments[:, 1]
    edge = points[ends1] - points[ends0]
    edge_length = np.linalg.norm(edge, axis=1)
    normal0 = normals[ends0]
    normal1 = normals[ends1]
    length0 = np.linalg.norm(normal0, axis=1)
    length1 = np.linalg.norm(normal1, axis=1)
    usable = (edge_length > 1e-12) & (length0 > 1e-12) & (length1 > 1e-12)
    direction = np.zeros_like(edge)
    direction[usable] = edge[usable] / edge_length[usable, None]
    along0 = np.abs(np.einsum("ij,ij->i", normal0, direction))
    along1 = np.abs(np.einsum("ij,ij->i", normal1, direction))
    along0[usable] /= length0[usable]
    along1[usable] /= length1[usable]
    mask[usable] = (along0[usable] >= alignment) & (along1[usable] >= alignment)
    return mask


TUNNEL_COLOR = (0.86, 0.15, 0.12)


def show_gabriel_graph(points, segments, tunnels):
    """Draw Gabriel edges. Tunnel edges, when marked, are a separate red network."""
    if len(segments) == 0:
        return
    if tunnels is None:
        ps.register_curve_network("Gabriel graph", points, segments, radius=0.0015)
        return
    ordinary = segments[~tunnels]
    marked = segments[tunnels]
    if len(ordinary):
        ps.register_curve_network("Gabriel graph", points, ordinary, radius=0.0015)
    if len(marked):
        network = ps.register_curve_network("tunnels", points, marked, radius=0.0015)
        network.set_color(TUNNEL_COLOR)


def main():
    positional, radius_factor, tunnel_alignment = parse_args(sys.argv[1:])
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
    tunnels = None
    edge_report = f"{graph.no_edges()} Gabriel edges"
    if tunnel_alignment is not None:
        tunnels = tunnel_edges(points, segments, tunnel_alignment)
        edge_report += (
            f", {int(np.count_nonzero(tunnels))} tunnel edges "
            f"at alignment {tunnel_alignment}"
        )
    print(
        f"{source}: {len(points)} points, radius factor {radius_factor}, "
        f"radius {radius:.6g}, {edge_report}"
    )

    ps.init()
    ps.set_ground_plane_mode("shadow_only")
    ps.register_point_cloud("points", points, radius=0.004)
    show_gabriel_graph(points, segments, tunnels)
    ps.show()


if __name__ == "__main__":
    main()
