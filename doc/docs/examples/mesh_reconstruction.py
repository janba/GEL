# No file arguments. Reconstructs a synthetic sphere and writes reconstructed.obj.
import pygel3d.hmesh as hmesh
import numpy as np

# Load or generate point cloud
# points: Nx3 numpy array
# normals: Nx3 numpy array

# Example usage (with synthetic data)
theta = np.linspace(0, 2*np.pi, 100)[:-1]
phi = np.linspace(0.1, np.pi-0.1, 50)
theta, phi = np.meshgrid(theta, phi)

x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

points = np.stack([x.flatten(), y.flatten(), z.flatten()], axis=1)
normals = points / np.linalg.norm(points, axis=1, keepdims=True)

m = hmesh.rsr_recon(points, normals, max_normal_ang=10)
hmesh.close_holes(m, max_size=1000)
hmesh.triangulate(m)
for _ in range(5):
    hmesh.maximize_min_angle(m)
hmesh.save("reconstructed.obj", m)
