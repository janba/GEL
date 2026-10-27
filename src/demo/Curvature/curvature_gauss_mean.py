"""Compare Gaussian and mean curvature on the torus.

The mesh is data/Solids/torus.obj: a torus of revolution about the Y axis
with major radius 1 and minor radius 1/2. Vertex normals point out of the
solid. Six scalar fields are shown in one window, one after another.
"""

from pygel3d import hmesh, gl_display as gl
from sys import argv
from numpy import array, hypot

# Major and minor radii of data/Solids/torus.obj.
R = 1.0
r = 0.5


def analytical_gauss_mean(m):
    """Analytical Gaussian and mean curvature at the vertex positions.

    With the outward normal, the principal curvatures are 1/r and
    (rho - R) / (r * rho), where rho is the distance to the Y axis.
    Gaussian curvature is their product and mean curvature is their average.
    """
    pos = m.positions()
    rho = hypot(pos[:, 0], pos[:, 2])
    k_parallel = (rho - R) / (r * rho)
    k_meridian = 1.0 / r
    return k_parallel * k_meridian, 0.5 * (k_parallel + k_meridian)


def show_scalar(viewer, mesh, title, values, reference=None):
    print(title)
    print("  range:", values.min(), values.max())
    if reference is not None:
        error = abs(values - reference)
        print("  mean abs error:", error.mean(), " max abs error:", error.max())
    viewer.display(mesh, mode='s', data=values, bg_col=(1, 1, 1))


fn = "../../../data/Solids/torus.obj" if len(argv) < 2 else argv[1]
m = hmesh.load(fn)
hmesh.triangulate(m)

K_ref, H_ref = analytical_gauss_mean(m)
K_bary = array([m.gaussian_curvature(v, use_mixed_area=False) for v in m.vertices()])
K_mixed = array([m.gaussian_curvature(v, use_mixed_area=True) for v in m.vertices()])
H_bary = array([m.mean_curvature(v, use_mixed_area=False) for v in m.vertices()])
H_mixed = array([m.mean_curvature(v, use_mixed_area=True) for v in m.vertices()])

viewer = gl.Viewer()
show_scalar(viewer, m, "Displaying analytical Gaussian curvature", K_ref)
show_scalar(viewer, m, "Displaying Gaussian curvature (barycentric area)", K_bary, K_ref)
show_scalar(viewer, m, "Displaying Gaussian curvature (mixed area)", K_mixed, K_ref)
show_scalar(viewer, m, "Displaying analytical mean curvature", H_ref)
show_scalar(viewer, m, "Displaying mean curvature (barycentric area)", H_bary, H_ref)
show_scalar(viewer, m, "Displaying mean curvature (mixed area)", H_mixed, H_ref)
