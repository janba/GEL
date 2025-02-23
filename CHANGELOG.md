Version 0.6.0
- Added type hints to Python interface
- Added Rotation System Reconstruction to GEL and PyGEL
- Vastly improved kD-Tree for m-nearest neighbor queries.
- Added complete mesh volume and area computation for HMesh
- updated C++ interface for HMesh. It was a bit random which functions were members of Manifold and which not. Now, all functions that operate on mesh elements are members. However, for backward compatibility the old functions are retained.