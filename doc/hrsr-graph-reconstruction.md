# Hierarchical RsR from the contraction graph

Two reconstructions are available and can be compared from PyGEL.

- `hrsr_recon` is the previous method. Collapse still builds a nearest-neighbor graph, contracts it, and keeps the surviving points and their averaged normals. The edges are discarded. RsR then builds a new nearest-neighbor graph and triangulates that.
- `hrsr_recon_graph` keeps the contracted graph and gives it to RsR. The normals on that graph are the contraction averages, normalized to unit length. In C++ the entry points are `point_cloud_collapse_reexpand` and `point_cloud_collapse_reexpand_graph`. The shared reconstruction is `graph_to_mesh`.

RsR still builds a spanning tree of the supplied graph, adds the remaining graph edges that pass the rotation-system checks, and then fills triangles from the tree. On the graph path that fill runs even when `genus` is 0, so no handles are requested. A new edge created by the fill may be as long as the nearest-neighbor radius (`num_neighbors`). The edges that already exist are only those of the contraction graph.

## Where this leaves the new reconstruction

The contraction graph is connected, but it does not contain enough edges to describe the surface. It is seeded with about five neighbors per point. After contraction it is only a little denser than its spanning tree. On `owl-little` (10k points, one collapse iteration) the tree has 4819 edges and the graph has 6981.

`hrsr_recon` covers the surface because the rebuilt graph has 30–70 neighbors per point (about 44k edges on the same cloud). `hrsr_recon_graph` cannot do that from the edges it receives. Filling the tree removes the severe fragmentation, and on that cloud the coarse mesh is one piece of about 4800 vertices plus a handful of tiny components. Reexpansion failures drop sharply, but the mesh is still poorer than the nearest-neighbor reconstruction. The limit is the missing connections in the simplified graph, not the spanning tree or the reexpansion.
