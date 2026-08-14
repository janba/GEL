# Usage: python load_process_save.py <input.obj> <output.obj>
import pygel3d.hmesh as hmesh
from sys import argv

def print_stats(m):
    # Print statistics
    print(f"Vertices: {len(m.vertices())}")
    print(f"Faces: {len(m.faces())}")
    print(f"Valid: {hmesh.valid(m)}")
    print(f"Closed: {hmesh.closed(m)}")

# Load mesh
m = hmesh.load(argv[1])
print_stats(m)

# Close holes and clean
hmesh.close_holes(m)
hmesh.triangulate(m)
hmesh.remove_caps(m)
hmesh.remove_needles(m)

# Print statistics after cleaning
print("After cleaning:")
print_stats(m)

# Smooth
hmesh.taubin_smooth(m,30)

# Save result
hmesh.save(argv[2], m)
