from pygel3d.hmesh import rsr_recon, Manifold
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

vertices, normals = obj_load('../../../data/PointClouds/owl-lines.obj')

m = rsr_recon(vertices,normals,False)
viewer = gl.Viewer()
viewer.display(m)