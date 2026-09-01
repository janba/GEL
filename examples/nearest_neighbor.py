# No file arguments. Builds a kD-tree on random points and times queries.
from pygel3d.spatial import I3DTree
import numpy as np
import time

# Generate random point cloud
n_points = 1000000
points = np.random.rand(n_points, 3) * 100

# Build kD-tree
print("Building tree...")
start = time.time()
tree = I3DTree()
for i, p in enumerate(points):
    tree.insert(p, i)
build_time = time.time() - start
print(f"Built tree with {n_points} points in {build_time:.2f}s")

# Query nearest neighbors
n_queries = 100000
query_points = np.random.rand(n_queries, 3) * 100

start = time.time()
for q in query_points:
    nearest_idx = tree.closest_point(q,1e20)
query_time = time.time() - start

print(f"Performed {n_queries} queries in {query_time:.2f}s")
print(f"Average query time: {query_time/n_queries*1000000:.2f}us")
