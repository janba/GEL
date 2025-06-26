# GEL Contributors
The following is the list of contributors. It is hopefully complete, but please let me know if you have contributed without being acknowldged or if you would like your specific contribution to be specified.

- Andreas Bærentzen created the library and has written most of the code not contributed by the people below.
- Jeppe Revall Frisvad has contributed greatly to the CGLA part of GEL and elsewhere. He particularly worked on the quaternion classes.
- Christian Thode Larsen rewrote the classes in HMesh, switching from lists to vectors for mesh attributes.
- Anders Wang Kristensen suggested important improvements to HMesh and wrote the GL Console code.
- Emil Gæde - added the multilevel LS skeletonization and graph components crucuial for Rotation System Reconstruction. This include tracking of connecting components in a graph when the graph is modified.
- Ruiqi Cui ported the Rotation System Reconstruction implementation to GEL.
- Karran Pandey added the original code for Skeleton to Face Extrusion Quad Meshes.
- Rasmus Paulsen, Bjarke Jakobsen, Marek Misztal helped with structure and compilation on various platforms with testing and ideas.
- Henrik Aanæs provided the original LinAlg package, which is a set of simple but very useful Lapack bindings. However, these were later removed. GEL and PyGEL do not require a numerics library, and we now use Eigen in applications that use GEL.
- Bent Dalgaard made early contributions: specifically the BSP Tree.
- Florian Gawrilowicz
- Morten Nobel-Jørgensen
- Kasper Steenstrup
- Jakob Wilm
- Tim Felle Olsen 
- Johnny Nunez
- Leonardo Mariga
- Tim Player
- rasmuschriste