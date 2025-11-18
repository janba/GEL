from sys import argv
from pygel3d.hmesh import rsr_recon, Manifold, save, flip_orientation
from pygel3d import gl_display as gl

def obj_load(file_path):
    vertices = []  # List to store vertex coordinates (x, y, z)
    normals = []  # List to store vertex normals (nx, ny, nz)

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if not parts:
                continue

            # Read vertex coordinates (v x y z)
            if parts[0] == 'v':
                vertices.append(tuple(map(float, parts[1:4])))

            # Read vertex normals (vn nx ny nz)
            elif parts[0] == 'vn':
                normals.append(tuple(map(float, parts[1:4])))

    return vertices, normals

viewer = gl.Viewer()

if len(argv) > 1:
    # Load from command line argument
    vertices, normals = obj_load(argv[1])
    m = rsr_recon(vertices,normals, use_Euclidean_distance=True)
    flip_orientation(m)
    viewer.display(m, smooth=False, mode='g')
    save("out.obj", m)

else:
    # larger point cloud.
    vertices, normals = obj_load('../../../data/PointClouds/owl-little.obj')
    m = rsr_recon(vertices,normals)
    flip_orientation(m)
    viewer.display(m, smooth=False, mode='g')
    save("owl.obj", m)

    # Object with non-zero genus
    vertices, normals = obj_load('../../../data/PointClouds/Capital_A.obj')
    m = rsr_recon(vertices, normals, True, k=30, n=40)
    viewer.display(m, smooth=False, mode='g', reset_view=True)
    save("A.obj", m)