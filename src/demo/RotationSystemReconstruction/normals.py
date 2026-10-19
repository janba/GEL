""" Show oriented normals for a point cloud.

    python normals.py [cloud.obj] [cluster|mst|spt|both|all]

    mst propagates signs along the tunnel-weighted spanning forest of the
    radius-limited Delaunay graph. rsr_recon uses that same alignment.
    spt propagates along a shortest-path tree of the k-nearest neighbor
    graph. cluster joins tangent edges
    before gap edges. both registers mst and cluster. all registers all
    three. The mst and spt views draw their trees. Edges that run the
    sheet cone check are red. Other MST edges are blue-gray, and other
    SPT edges are green.
    The default cloud is owl-lines and the default orientation is mst.
"""
from sys import argv

import numpy as np
import polyscope as ps
from pygel3d.hmesh import rsr_normals, rsr_orientation_tree


def load_points(path):
    points = []
    with open(path) as cloud:
        for line in cloud:
            parts = line.split()
            if parts and parts[0] == 'v':
                points.append(tuple(map(float, parts[1:4])))
    return np.asarray(points, dtype=np.float64)


def parse_args(args):
    methods = {'cluster', 'mst', 'spt', 'both', 'all'}
    path = '../../../data/PointClouds/owl-lines.obj'
    orientation = 'mst'
    for arg in args:
        if arg in methods:
            orientation = arg
        else:
            path = arg
    return path, orientation


path, orientation = parse_args(argv[1:])
points = load_points(path)
# rsr.py reconstructs this directory's clouds with Euclidean distance.
print(f'{len(points)} points from {path}, orientation {orientation}')

TREE_COLOR = {
    'mst': (0.45, 0.58, 0.72),
    'spt': (0.18, 0.62, 0.48),
}
CONE_COLOR = (0.93, 0.16, 0.12)


def show_orientation_tree(points, edges, cone, method):
    if len(edges) == 0:
        print(f'orientation {method} has no edges')
        return
    span = np.linalg.norm(points.max(axis=0) - points.min(axis=0))
    radius = 0.0015 * span if span > 0 else 0.001
    network = ps.register_curve_network(f'orientation {method}', points, edges, radius=radius)
    colors = np.tile(np.array(TREE_COLOR[method], dtype=np.float64), (len(edges), 1))
    colors[np.asarray(cone, dtype=bool)] = CONE_COLOR
    network.add_color_quantity("cone check", colors, defined_on='edges', enabled=True)
    print(f'{len(edges)} {method.upper()} edges, {int(np.count_nonzero(cone))} run the cone check')


def add_oriented(cloud, points, method, enabled):
    if method in TREE_COLOR:
        normals, edges, cone = rsr_orientation_tree(
            points, use_Euclid_dist=True, orientation=method)
        cloud.add_vector_quantity(method, normals, enabled=enabled)
        show_orientation_tree(points, edges, cone, method)
        return
    cloud.add_vector_quantity(
        method,
        rsr_normals(points, use_Euclid_dist=True, orientation=method),
        enabled=enabled)


ps.init()
ps.set_ground_plane_mode("shadow_only")
cloud = ps.register_point_cloud("points", points)
if orientation == 'all':
    add_oriented(cloud, points, 'mst', True)
    add_oriented(cloud, points, 'spt', False)
    add_oriented(cloud, points, 'cluster', False)
elif orientation == 'both':
    add_oriented(cloud, points, 'mst', True)
    add_oriented(cloud, points, 'cluster', False)
else:
    add_oriented(cloud, points, orientation, True)
ps.show()
