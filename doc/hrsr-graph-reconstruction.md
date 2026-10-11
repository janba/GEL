# Hierarchical RsR from the contraction graph

Two reconstructions are available and can be compared from PyGEL.

- `hrsr_recon` collapses a 5-neighbor graph, discards the edges, and builds a new nearest-neighbor graph from the collapsed points.
- `hrsr_recon_graph` builds the collapse graph the way RsR would, then gives the surviving graph to RsR. In C++ the entry points are `point_cloud_collapse_reexpand` and `point_cloud_collapse_reexpand_graph`. The shared reconstruction is `graph_to_mesh`.

In Python the functions are `hmesh.hrsr_recon` and `hmesh.hrsr_recon_graph`. The arguments match `rsr_recon`, with `collapse_iters` (default 1) and `skip_reexpansion` (default false).

## Normals and orientation

Both paths share the normal estimate. Normals passed in by the caller are kept and normalized. An empty normal array is filled in before the collapse.

At each point the estimator scores neighborhoods of 8, 12, 18, 27, 40, 60, 90, 128, and 192 neighbors, taken as prefixes of one distance-sorted search. It keeps the planar covariance with the lowest eigenentropy. A prefix counts as planar when its middle eigenvalue is at least 0.15 times the largest. The two in-plane eigenvalues are averaged before the entropy, so a long thin neighborhood is not preferred over a round one. Scores within `1e-3` keep the smaller neighborhood. If no prefix is planar, the most planar prefix is kept. The sign of the normal is arbitrary. The run prints `adaptive normal median k`. A smooth sparse region stays on a few neighbors. A noisy region grows until the noise averages out.

Orientation propagates along a minimum spanning tree of that neighborhood, built on the smoothed positions. The edge cost is `1 - |n·m|` plus `(|e·n| + |e·m|) / scale`, where `scale` is the point's farthest-neighbor distance. The second term is omitted when `scale` is below `1e-12`. An edge that leaves the tangent plane is expensive, which keeps the tree from tunneling through the cloud. A vertex is flipped when its normal dots negatively with its parent. The smoothed positions are then discarded. Collapse and reconstruction use the original points.

## Collapse graph

The 5-neighbor graph exists so that collapse only merges close points. Contraction rewires those edges and does not search for new ones, so after one collapse on `owl-little` it had 6981 edges against a spanning tree of 4819. That is enough to stay connected and not enough for RsR, which can only turn an edge it is given into a face.

`hrsr_recon_graph` therefore seeds the full neighborhood (`num_neighbors`, with the same normal-angle test as RsR). An edge of that graph may be contracted only when its Euclidean length is within the one-ring of an endpoint. The one-ring is the distance to the `initial_neighbors`-th nearest live point, 5 by default. On the first iteration that radius is the distance of that neighbor in the seed search. Later iterations search again, because the points have moved. Edges outside the ring stay in the graph. They live in a side list and are absent from the contraction queue, which holds the one-ring. Export still gives RsR every edge that was not contracted.

The collapse queue is ordered by tangent distance times the sum of the vertex weights, including when reconstruction itself uses Euclidean distance. Tangent distance is not a length cap: an edge that jumps along the normals can look short in that measure, so the one-ring test is what keeps it out of the queue. `hrsr_recon` uses the same queue key and still collapses the 5-neighbor seed.

On export, RsR receives the survivors, the averaged normals normalized to unit length, and every edge that was not contracted. Face edges still have to be shorter than the neighbor at rank `num_neighbors * 2/3`, the same limit as a point-cloud reconstruction. The spanning tree and the face order use that length. The spanning tree is not filled when `genus` is 0.

A stored neighbor list is reused only when it is an exact search of the requested count. The kD-tree heap drops distance ties once it is full, so a prefix of a larger search is a different neighbor set. A rank distance, such as the one-ring radius or the `2/3` length cap, may be read from a longer search.

## Owl-little, one collapse, Euclidean, genus 0, reexpansion skipped

The coarse meshes are the same size. The times were measured while every reconstruction edge still sat in the contraction star: about 1.7 s on the graph path against 0.15 s for this cloud of 10k points. The queue now holds only the one-ring. That run has not been repeated. The mesh counts below are from the skipped-reexpansion comparison.

| | vertices | faces | boundaries | graph edges | face candidates |
| --- | --- | --- | --- | --- | --- |
| `hrsr_recon` | 5001 plus one triangle | 9754 | 11 | rebuilt nearest-neighbor graph | 101032 |
| `hrsr_recon_graph` | 5004 | 9754 | 10 | 105312 | 93840 |

## Reexpansion

Reexpansion is on by default. `skip_reexpansion=True` leaves the coarse mesh, which is the mesh in the table above. On the clouds tried so far that coarse mesh is the better surface. Reinsertion is what damages it.

Each coarse vertex sits at `v_bar`, the weighted average stored for that collapse. Reexpansion replays the log in reverse and finds the vertex by exact `double` equality with `v_bar`. For each record it does four things:

1. It flips link edges that the segment between the two original samples passes near, up to 100 times.
2. It moves the coarse vertex from `v_bar` to the active sample.
3. It splits that vertex. The original vertex stays at the active sample. The new vertex is created at the same place and is then moved to the latent sample.
4. It flips edges around the new vertex, longest first, for `refinement_iterations` (default 1). A corner under 22.5° is flipped when the four corners beside the edge are acute and the flip is legal. The flip does not have to improve the triangles.

The split gap minimizes a dihedral cost that is zero on a flat surface, plus a minimum-angle cost. The angle cost is the smallest corner's deviation from 60°, times the shortest edge. A corner under 3° adds `min_angle_threshold_penalty` (0.1) on top of that deviation. The penalty used to replace the deviation, and because 0.1 is smaller than the deviation of an ordinary sharp corner, those splits preferred slivers. The score is still multiplied by the shortest edge, so a sliver whose short edge is tiny can remain cheap.

Orientation repairs for this phase are present in `reexpand_points` and compiled out (`if (false)`). They delete boundary edges and flip spokes whose triangles face against the vertex normal. Enabling them removed edges that the refinement sort still held, and reexpansion aborted on `manifold.in_use`.

A crossing repair that cannot flip leaves the coarse vertex at the active sample and skips the split. The position map is updated only after a successful split, so the rest of that vertex's collapse chain is then missed.
