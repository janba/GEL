# Hierarchical RsR from the contraction graph

Two reconstructions are available and can be compared from PyGEL.

- `hrsr_recon` collapses a 5-neighbor graph, discards the edges, and builds a new nearest-neighbor graph from the collapsed points.
- `hrsr_recon_graph` builds the collapse graph the way RsR would, then gives the surviving graph to RsR. In C++ the entry points are `point_cloud_collapse_reexpand` and `point_cloud_collapse_reexpand_graph`. The shared reconstruction is `graph_to_mesh`.

Both paths share the normal estimate. Normals passed in by the caller are kept and normalized. An empty normal array is filled in before the collapse. At each point the estimator scores neighborhoods of 8, 12, 18, 27, 40, 60, 90, 128, and 192 neighbors. It keeps the planar covariance with the lowest eigenentropy, and the smaller neighborhood when those scores tie. The two in-plane eigenvalues are averaged before the entropy, so a long thin neighborhood is not preferred over a round one. A smooth sparse region stays on a few neighbors. A noisy region grows until the noise averages out. Orientation still propagates along the minimum spanning tree of the neighborhood, and an edge that leaves the tangent plane is expensive.

The 5-neighbor graph exists so that collapse only merges close points. Contraction rewires those edges and does not search for new ones, so after one collapse on `owl-little` it had 6981 edges against a spanning tree of 4819. That is enough to stay connected and not enough for RsR, which can only turn an edge it is given into a face.

`hrsr_recon_graph` therefore seeds the full neighborhood (`num_neighbors`, with the same normal-angle test as RsR). An edge of that graph may be contracted only when its Euclidean length is within the one-ring of an endpoint. The one-ring is the distance to the `initial_neighbors`-th nearest live point, 5 by default, and it is recomputed at each collapse iteration. Edges outside that ring stay in the graph.

The collapse queue is ordered by tangent distance times the sum of the vertex weights. Tangent distance is not a length cap: an edge that jumps along the normals can look short in that measure, so the one-ring test is what keeps it out of the queue. `hrsr_recon` uses the same queue key and still collapses the 5-neighbor seed.

On export, RsR receives the survivors, the averaged normals normalized to unit length, and every edge that was not contracted. Face edges still have to be shorter than the neighbor at rank `num_neighbors * 2/3`, the same limit as a point-cloud reconstruction. The spanning tree and the face order use that length. The spanning tree is not filled when `genus` is 0.

## Owl-little, one collapse, Euclidean, genus 0, reexpansion skipped

The coarse meshes are the same size. Collapse on the graph path is slower because each merge rewrites a much larger star: about 1.7 s against 0.15 s for this cloud of 10k points.

| | vertices | faces | boundaries | graph edges | face candidates |
| --- | --- | --- | --- | --- | --- |
| `hrsr_recon` | 5001 plus one triangle | 9754 | 11 | rebuilt nearest-neighbor graph | 101032 |
| `hrsr_recon_graph` | 5004 | 9754 | 10 | 105312 | 93840 |
