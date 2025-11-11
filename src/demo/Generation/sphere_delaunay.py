from pygel3d import hmesh, gl_display as gl 
import numpy as np

def generate_random_sphere_points(n=64):
    """
    Generate n uniformly distributed random points on the unit sphere.
    Uses the method of generating random points in 3D and normalizing.
    """
    # Generate random points from a normal distribution
    points = np.random.randn(n, 3)
    
    # Normalize each point to lie on the unit sphere
    norms = np.linalg.norm(points, axis=1, keepdims=True)
    points = points / norms
    
    return points

# Generate 64 random points on the unit sphere
points = generate_random_sphere_points(64)
m = hmesh.sphere_delaunay(points)

# Perturb vertices
pos = m.positions()
for v in m.vertices():
    pos[v] *= np.random.uniform(0.5, 1.5)

# View noisy mesh
v = gl.Viewer()
v.display(m) # HIT ESCAPE TO CONTINUE

# Smooth and view again
hmesh.anisotropic_smooth(m, no_iters=10)
v.display(m) # HIT ESCAPE TO EXIT
