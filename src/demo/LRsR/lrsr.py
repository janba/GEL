"""Local rotation system reconstruction.

    python lrsr.py [--radius FACTOR] [points.obj] [out.obj]

With no file, 400 points on the unit sphere are reconstructed and shown.
The window draws the mesh and the tunnel-weighted spanning forest of the
same Delaunay graph. Modes 1 and 3 build that forest and use its edges.
Mode 2 does not. `kGrowthMode` in LRsR.cpp selects the mode and is 2.
The radius is FACTOR times the median nearest-neighbor distance. The default
factor is 4. An .obj point cloud uses each line whose first word is v.
The mesh is written to out.obj when no output path is given.
Reconstruction prints `LRsR mode: M` and `LRsR iterations: N`. N is how
many triangles the queue accepted. Mode 2 does not count its seed triangle.
Mode 1 queues tree edges and edges of accepted triangles together. Mode 2
starts from the lowest-score triangle and spreads along a front. Mode 3
exhausts the tree edges first, then spreads from the triangles already
added. A queued triangle is scored by its projected perimeter times
`2 - cos(alpha)` for each face it already touches. A popped entry whose
score has risen is put back at that score. Mode 2 also winds each queued
triangle against the face across the edge that queued it, and drops a
popped triangle that disagrees with a neighbor. Each connected component is
then oriented from its first triangle, using only the cyclic order of
the vertices.
"""

import sys
from collections import defaultdict
from os.path import splitext

import numpy as np
from pygel3d import gl_display as gl
from pygel3d.graph import delaunay_edges, minimum_spanning_forest
from pygel3d.hmesh import count_boundary_curves, lrsr_recon, save, valid


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


def sanity(points, mesh, radius_factor):
    """Forest edges absent from the mesh, non-Delaunay faces, and bad fans.

    A fan is fine when the triangles at one input point are a single cycle
    or a set of paths. A branch, or a cycle plus anything else, is not.
    The returned forest is the tunnel-weighted spanning forest of the same
    Delaunay graph. It is drawn beside the mesh and is not used to seed the
    fill. An empty mesh has no positions array. Every forest edge is then
    missing, and every input point is without a face.
    """
    if len(mesh.vertices()) == 0:
        radius = radius_factor * median_neighbor_spacing(points)
        graph = delaunay_edges(points, radius)
        forest, *_ = minimum_spanning_forest(graph, True)
        forest_edges = 0
        for index in range(len(points)):
            for neighbor in forest.neighbors(index):
                if index < int(neighbor):
                    forest_edges += 1
        return {
            "missing_forest": forest_edges,
            "forest_edges": forest_edges,
            "non_delaunay_faces": 0,
            "bad_fans": 0,
            "closed": 0,
            "open_fans": 0,
            "points_without_faces": len(points),
            "forest": forest,
        }
    positions = np.array(mesh.positions())

    def key_of(point):
        return (round(float(point[0]), 8), round(float(point[1]), 8), round(float(point[2]), 8))

    point_at = defaultdict(list)
    for index, point in enumerate(points):
        point_at[key_of(point)].append(index)
    point_id = {}
    for vertex in mesh.vertices():
        hits = point_at.get(key_of(positions[int(vertex)]), [])
        if len(hits) == 1:
            point_id[int(vertex)] = hits[0]

    radius = radius_factor * median_neighbor_spacing(points)
    graph = delaunay_edges(points, radius)
    delaunay = set()
    for index in range(len(points)):
        for neighbor in graph.neighbors(index):
            delaunay.add(tuple(sorted((int(index), int(neighbor)))))
    forest, *_ = minimum_spanning_forest(graph, True)
    missing_forest = 0
    forest_edges = 0
    mesh_edges = set()
    non_delaunay_faces = 0
    incident = defaultdict(list)
    for face in mesh.faces():
        vertices = [int(vertex) for vertex in mesh.circulate_face(face, "v")]
        if any(vertex not in point_id for vertex in vertices):
            continue
        ids = [point_id[vertex] for vertex in vertices]
        bad = False
        for start, end in ((0, 1), (1, 2), (2, 0)):
            edge = tuple(sorted((ids[start], ids[end])))
            mesh_edges.add(edge)
            if edge not in delaunay:
                bad = True
        if bad:
            non_delaunay_faces += 1
        for corner in range(3):
            incident[ids[corner]].append((ids[(corner + 1) % 3], ids[(corner + 2) % 3]))
    for index in range(len(points)):
        for neighbor in forest.neighbors(index):
            if index < int(neighbor):
                forest_edges += 1
                if (index, int(neighbor)) not in mesh_edges:
                    missing_forest += 1

    bad_fans = 0
    closed = 0
    open_fans = 0
    for triangles in incident.values():
        adjacency = defaultdict(set)
        for left, right in triangles:
            adjacency[left].add(right)
            adjacency[right].add(left)
        seen = set()
        components = 0
        branched = False
        cyclic = False
        for start in list(adjacency):
            if start in seen:
                continue
            components += 1
            stack = [start]
            seen.add(start)
            nodes = []
            edge_ends = 0
            while stack:
                node = stack.pop()
                nodes.append(node)
                if len(adjacency[node]) > 2:
                    branched = True
                for other in adjacency[node]:
                    edge_ends += 1
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
            if edge_ends // 2 >= len(nodes) and len(nodes) >= 3:
                cyclic = True
        if branched or (cyclic and components > 1):
            bad_fans += 1
        elif cyclic and components == 1 and not branched:
            closed += 1
        else:
            open_fans += 1
    return {
        "missing_forest": missing_forest,
        "forest_edges": forest_edges,
        "non_delaunay_faces": non_delaunay_faces,
        "bad_fans": bad_fans,
        "closed": closed,
        "open_fans": open_fans,
        "points_without_faces": len(points) - len(incident),
        "forest": forest,
    }


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
    output = positional[1] if len(positional) > 1 else "out.obj"

    mesh = lrsr_recon(points, radius_factor)
    n_vertices = len(mesh.vertices())
    n_faces = len(mesh.faces())
    check = sanity(points, mesh, radius_factor)
    print(
        f"{source}: {len(points)} points, radius factor {radius_factor}, "
        f"{n_vertices} vertices, {n_faces} faces, "
        f"{count_boundary_curves(mesh)} boundary curves, "
        f"valid={valid(mesh)}"
    )
    print(
        f"sanity: missing forest edges {check['missing_forest']}"
        f"/{check['forest_edges']}, "
        f"non-Delaunay faces {check['non_delaunay_faces']}, "
        f"bad fans {check['bad_fans']}, "
        f"closed umbrellas {check['closed']}, "
        f"open fans {check['open_fans']}, "
        f"points without faces {check['points_without_faces']}"
    )
    save(output, mesh)
    gl.Viewer().display(mesh, check["forest"], smooth=False, mode="n")


if __name__ == "__main__":
    main()
