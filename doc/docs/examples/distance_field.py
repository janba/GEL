# Expects bunny.obj in the current directory.
# Writes distance_field.npy.
from sys import argv
from pygel3d import hmesh
import numpy as np

# Load mesh
m = hmesh.load(argv[1])

# Create distance object
dist = hmesh.MeshDistance(m)

# Create sampling grid
grid_res = 50
bbox_min, bbox_max = hmesh.bbox(m)

# Expand bbox slightly
margin = 0.1
bbox_min = [x - margin for x in bbox_min]
bbox_max = [x + margin for x in bbox_max]

# Sample points
x = np.linspace(bbox_min[0], bbox_max[0], grid_res)
y = np.linspace(bbox_min[1], bbox_max[1], grid_res)
z = np.linspace(bbox_min[2], bbox_max[2], grid_res)

# Compute distance field
distances = np.zeros((grid_res, grid_res, grid_res))
for i, xi in enumerate(x):
    for j, yj in enumerate(y):
        for k, zk in enumerate(z):
            point = [xi, yj, zk]
            distances[i, j, k] = dist.signed_distance(point)

print(f"Distance field computed: {distances.shape}")
print(f"Min distance: {distances.min():.4f}")
print(f"Max distance: {distances.max():.4f}")

# Save distance field (e.g., as numpy array)
np.save("distance_field.npy", distances)
