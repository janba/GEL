"""Local rotation system reconstruction.

    python lrsr.py [--radius FACTOR] [points.obj] [out.obj]

With no file, 400 points on the unit sphere are reconstructed and shown.
The radius is FACTOR times the median nearest-neighbor distance. The default
factor is 4. An .obj point cloud uses each line whose first word is v.
The mesh is written to out.obj when no output path is given.
"""

import sys
from os.path import splitext

import numpy as np
from pygel3d import gl_display as gl
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
    print(
        f"{source}: {len(points)} points, radius factor {radius_factor}, "
        f"{n_vertices} vertices, {n_faces} faces, "
        f"{count_boundary_curves(mesh)} boundary curves, "
        f"valid={valid(mesh)}"
    )
    save(output, mesh)
    gl.Viewer().display(mesh, smooth=False, mode="g")


if __name__ == "__main__":
    main()
