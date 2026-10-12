""" Show the normals that rotation-system reconstruction computes for a point cloud. """
from sys import argv

import numpy as np
import polyscope as ps
from pygel3d.hmesh import rsr_normals


def load_points(path):
    points = []
    with open(path) as cloud:
        for line in cloud:
            parts = line.split()
            if parts and parts[0] == 'v':
                points.append(tuple(map(float, parts[1:4])))
    return np.asarray(points, dtype=np.float64)


path = argv[1] if len(argv) > 1 else '../../../data/PointClouds/owl-lines.obj'
points = load_points(path)
# rsr.py reconstructs this directory's clouds with Euclidean distance.
normals = rsr_normals(points, use_Euclid_dist=True)

ps.init()
ps.set_ground_plane_mode("shadow_only")
cloud = ps.register_point_cloud("points", points)
cloud.add_vector_quantity("normals", normals, enabled=True)
ps.show()
