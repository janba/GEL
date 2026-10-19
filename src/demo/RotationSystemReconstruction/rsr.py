from sys import argv, exit
from pygel3d.hmesh import rsr_recon, hrsr_recon, hrsr_recon_graph, save, flip_orientation, load, close_holes, triangulate
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

def ground_truth_cloud(file_path):
    mesh = load(file_path)
    if mesh is None:
        print(f"Could not load {file_path}")
        return None
    ids = list(mesh.vertices())
    if len(ids) == 0:
        print(f"{file_path} has no vertices")
        return None
    positions = mesh.positions()
    vertices = []
    normals = []
    for vid in ids:
        vertices.append(tuple(positions[vid]))
        normals.append(tuple(mesh.vertex_normal(vid)))
    return vertices, normals

def flag_name(arg):
    return arg.lower().lstrip('-')

viewer = gl.Viewer()

if len(argv) > 1:
    args = argv[1:]
    hierarchical = any(flag_name(arg) == 'h' for arg in args)
    gt_normals = any(flag_name(arg) == 'gt-normals' for arg in args)
    files = [arg for arg in args if flag_name(arg) not in ('h', 'gt-normals')]
    if not files:
        print("Missing input file")
        exit(1)
    filename = files[0]
    if gt_normals:
        loaded = ground_truth_cloud(filename)
        if loaded is None:
            exit(1)
        vertices, normals = loaded
        print(f"Loaded {len(vertices)} vertices and computed {len(normals)} vertex normals from {filename}")
    else:
        vertices, normals = obj_load(filename)
        print(f"Loaded {len(vertices)} vertices and {len(normals)} normals from {filename}")

    m = None
    if hierarchical:
        m = hrsr_recon_graph(vertices, normals, num_neighbors=50, max_neighbor_dist = 5, use_Euclid_dist=True) 
    else:
        m = rsr_recon(vertices, normals,  max_neighbor_dist=10, use_Euclid_dist=True, genus=1)
    close_holes(m, max_size=10000)
    triangulate(m, clip_ear=False)
    print("Reconstruction completed.")
    viewer.display(m, smooth=False, mode='g')
    save("out.obj", m)

else:
    # larger point cloud.
    vertices, normals = obj_load('../../../data/PointClouds/owl-lines.obj')
    m = rsr_recon(vertices, normals, use_Euclid_dist=True)
    # m = hrsr_recon_graph(vertices, normals, use_Euclid_dist=True, genus=0, collapse_iters=3, num_neighbors=40)
    flip_orientation(m)
    # viewer.display(m, smooth=False, mode='g')
    save("owl.obj", m)

    # # Object with non-zero genus
    # vertices, normals = obj_load('../../../data/PointClouds/Capital_A.obj')
    # m = hrsr_recon_graph(vertices, normals, use_Euclid_dist=True, genus=1, max_handle_dist=10)
    # viewer.display(m, smooth=False, mode='g', reset_view=True)
    # save("A.obj", m)