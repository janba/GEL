from sys import argv
from pygel3d.hmesh import rsr_recon, hrsr_recon, save, flip_orientation
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
    # Check if we should use the hierarchical method
    hierarchical = False
    if len(argv) > 2 and ( argv[1].lower() == '-h' or argv[2].lower() == '-h' ):
        hierarchical = True

    filename = argv[1] if not hierarchical else argv[2]
    # Load from command line argument
    vertices, normals = obj_load(filename)
    print(f"Loaded {len(vertices)} vertices and {len(normals)} normals from {filename}")

    m = None
    if hierarchical:
        m = hrsr_recon(vertices, normals)
    else:
        m = rsr_recon(vertices, normals)

    print("Reconstruction completed.")
    viewer.display(m, smooth=False, mode='g')
    save("out.obj", m)

else:
    # larger point cloud.
    vertices, normals = obj_load('../../../data/PointClouds/owl-little.obj')
    m = hrsr_recon(vertices, normals, use_Euclid_dist=True, genus=0)
    flip_orientation(m)
    viewer.display(m, smooth=False, mode='g')
    save("owl.obj", m)

    # Object with non-zero genus
    vertices, normals = obj_load('../../../data/PointClouds/Capital_A.obj')
    m = rsr_recon(vertices, normals, use_Euclid_dist=True, genus=1, max_handle_dist=10)
    viewer.display(m, smooth=False, mode='g', reset_view=True)
    save("A.obj", m)