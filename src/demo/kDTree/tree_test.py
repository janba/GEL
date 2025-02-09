from pygel3d import spatial
from random import uniform
import time

for trial in range(7):
    T = spatial.I3DTree()

    print("Building GEL I3DTree")
    pts = [(uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)) for _ in range(1000000)]
    for i,p in enumerate(pts):
        T.insert(p,i)
    T.build()

    if trial == 0:
        print("Validating I3DTree")
        for i,p in enumerate(pts):
            k, v = T.closest_point(p,1e-6)
            assert v == i

    print("Timing queries")
    start_time = time.time()
    indices = []
    for p in [(uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)) for _ in range(100000)]: 
        k, v = T.closest_point(p,1e32)
        indices.append(k)
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")