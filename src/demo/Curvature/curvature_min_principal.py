"""Compare minimum principal curvature on the torus.

The mesh is data/Solids/torus.obj: a torus of revolution about the Y axis
with major radius 1 and minor radius 1/2. Vertex normals point out of the
solid. Four scalar fields are shown in one window, one after another:
the analytical value, the one-ring paraboloid fit, the same fit with the
opposite vertex of each link edge, and that fit with Gaussian weights.
"""

from pygel3d import hmesh, gl_display as gl
from sys import argv
from numpy import array, hypot

# Major and minor radii of data/Solids/torus.obj.
R = 1.0
r = 0.5


def analytical_min_curvature(m):
    """Analytical minimum principal curvature at the vertex positions.

    With the outward normal this is the parallel curvature
    (rho - R) / (r * rho), where rho is the distance to the Y axis.
    The meridional curvature 1/r is never smaller.
    """
    pos = m.positions()
    rho = hypot(pos[:, 0], pos[:, 2])
    return (rho - R) / (r * rho)


def min_principal(m, use_opposite_vertices, use_weights):
    return array([
        m.principal_curvatures(v, use_opposite_vertices, use_weights)[0]
        for v in m.vertices()
    ])


def show_scalar(viewer, mesh, title, values, reference=None):
    print(title)
    print("  range:", values.min(), values.max())
    if reference is not None:
        error = abs(values - reference)
        print("  mean abs error:", error.mean(), " max abs error:", error.max())
    viewer.display(mesh, mode='s', data=values, bg_col=(1,1,1))


fn = "../../../data/Solids/torus.obj" if len(argv) < 2 else argv[1]
m = hmesh.load(fn)
hmesh.triangulate(m)

k_ref = analytical_min_curvature(m)
k_ring = min_principal(m, use_opposite_vertices=False, use_weights=False)
k_flaps = min_principal(m, use_opposite_vertices=True, use_weights=False)
k_weighted = min_principal(m, use_opposite_vertices=True, use_weights=True)

viewer = gl.Viewer()
show_scalar(viewer, m, "Displaying analytical minimum principal curvature", k_ref)
show_scalar(viewer, m, "Displaying minimum principal curvature (one ring)", k_ring, k_ref)
show_scalar(viewer, m, "Displaying minimum principal curvature (one ring + flaps)", k_flaps, k_ref)
show_scalar(viewer, m, "Displaying minimum principal curvature (one ring + flaps + weights)", k_weighted, k_ref)
