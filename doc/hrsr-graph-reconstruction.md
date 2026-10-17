# Hierarchical RsR from the contraction graph

Two reconstructions are available and can be compared from PyGEL.

- `hrsr_recon` collapses a 5-neighbor graph, discards the edges, and builds a new nearest-neighbor graph from the collapsed points.
- `hrsr_recon_graph` builds the collapse graph the way RsR would, then gives the surviving graph to RsR. In C++ the entry points are `point_cloud_collapse_reexpand` and `point_cloud_collapse_reexpand_graph`. The shared reconstruction is `graph_to_mesh`.

In Python the functions are `hmesh.hrsr_recon` and `hmesh.hrsr_recon_graph`. The arguments match `rsr_recon`, with `collapse_iters` (default 1) and `skip_reexpansion` (default false).

## Normals and orientation

Both paths share the normal estimate. Normals passed in by the caller are kept and normalized. An empty normal array is estimated and then oriented before the collapse.

At each point the estimator scores neighborhoods of 8, 10, 12, 14, 16, 20, 25, 30, 40, 50, 64, 128, and 192 neighbors, taken as prefixes of one distance-sorted search. It keeps the planar covariance with the lowest score. A prefix counts as planar when its middle eigenvalue is at least 0.15 times the largest. The score is `H(m, m, λ3) + (ln 2 − H(λ1, λ2))`, where `m = (λ1 + λ2) / 2` and `H` is the Shannon entropy of a spectrum in nats. A zero eigenvalue contributes nothing. A round flat patch scores `ln 2`. A flat strip scores higher by the entropy missing from the in-plane pair, so a larger rounder prefix replaces it while that prefix stays a surface. Scores within `1e-3` keep the smaller neighborhood. If no prefix is planar, the most planar prefix is kept. The sign of the normal is arbitrary. The run prints `adaptive normal median k`. A smooth round region stays on a few neighbors. A noisy or line-like region grows until the patch rounds out or the noise averages out.

Reconstruction orients on the kNN graph of `num_neighbors`. Euclidean orientation uses the original points. Tangent orientation uses the smoothed positions, which are then discarded. Collapse and reconstruction use the original points. `hmesh.rsr_normals` returns this empty-normal path. `orientation="mst"` is the default and is what reconstruction uses. `orientation="spt"` propagates along a shortest-path tree of that same graph. `orientation="cluster"` joins tangent edges before gap edges. `src/demo/RotationSystemReconstruction/normals.py` draws one of the three, `both` (mst and cluster), or `all`. It draws the MST and the shortest-path tree.

An undirected edge is scored once. Let `e` be the unit edge, `a = |n·e|`, and `b = |m·e|`. Tangent evidence is `(1-a)(1-b)`. Tunnel evidence is `a*b`. The score is `|n·m|` times the larger of those two, so parallel and anti-parallel normals score the same. A tangent edge and a tunnel edge both score 1. A 45° edge scores less. A zero-length edge scores 0. The function is `orientation_cluster_evidence`. Its `signed_evidence` is positive when a tangent edge's normals already agree, and when a tunnel edge's normals already face apart.

Clustering joins the tangent edges first, highest `tangent * |n·m|` first. An edge is tangent when its tangent evidence is at least its tunnel evidence and that rank is above `1e-4`. Every join flips the smaller side so the two normals agree. The remaining edges are joined afterwards, cheapest `orientation_edge_weight` first, and those joins agree as well. A gap between two parts that already meet along the surface does not get a vote, so it cannot reverse a patch. The run prints `Orientation clustering done`.

The spanning tree uses the same kNN graph. Its edge cost is `orientation_edge_weight`, `1 - |n·m| + (|e·n| + |e·m|) / scale`, where `scale` is that point's farthest-neighbor distance and `e` is the edge vector. The second term is omitted when `scale` is below `1e-12`. It is the penalty for an edge that leaves the tangent plane, so a tunnel loses to a tangent edge and stays out of the tree unless it is the link between two otherwise separate parts. A child agrees with its parent. `RSROpts::cone_check` (default false, `cone_check=True` in Python) changes that on a tough edge, where both `|n·ê|` and `|m·ê|` are at least 1/2. The cone then opens from the child along its normal toward the parent, with a 45° half-angle, and stops at twice the edge length. Another sheet in that cone gives the child the opposite sign from the parent. The parent is not flipped. An empty cone keeps the agree-flip. Each connected component keeps its own sign, and each tree is rooted at the component's lowest index. The run prints `correct normal orientation (mst, cone)` when the check is on, and `Orientation MST done`.

The shortest-path tree is the other propagation, `orientation="spt"`. It builds that same graph and pays `orientation_edge_weight`, but the queue key is the path sum from the root. The child flip is the same, including the cone check. Reconstruction does not use this tree. The run prints `correct normal orientation (spt)` and `Orientation SPT done`.

## Collapse graph

The 5-neighbor graph exists so that collapse only merges close points. Contraction rewires those edges and does not search for new ones, so after one collapse on `owl-little` it had 6981 edges against a spanning tree of 4819. That is enough to stay connected and not enough for RsR, which can only turn an edge it is given into a face.

`hrsr_recon_graph` therefore seeds the full neighborhood (`num_neighbors`, with the same normal-angle test as RsR). An edge of that graph may be contracted only when its Euclidean length is within the one-ring of an endpoint. The one-ring is the distance to the `initial_neighbors`-th nearest live point, 5 by default. On the first iteration that radius is the distance of that neighbor in the seed search. Later iterations search again, because the points have moved. Edges outside the ring stay in the graph. They live in a side list and are absent from the contraction queue, which holds the one-ring. Export still gives RsR every edge that was not contracted.

The collapse queue is ordered by tangent distance times the sum of the vertex weights, including when reconstruction itself uses Euclidean distance. Tangent distance is not a length cap: an edge that jumps along the normals can look short in that measure, so the one-ring test is what keeps it out of the queue. `hrsr_recon` uses the same queue key and still collapses the 5-neighbor seed.

On export, RsR receives the survivors, the averaged normals normalized to unit length, and every edge that was not contracted. The first face pass still requires an edge to be shorter than the neighbor at rank `num_neighbors * 2/3` at both endpoints, the same conservative limit as a point-cloud reconstruction. That pass runs before handle insertion. After handles, a second shortest-first pass tries the remaining graph edges whose length is still within each endpoint's longest angle-accepted neighbor. The same pass runs when `genus` is 0. Chord filling runs for every genus, including 0. Its chords join successive neighbors in the rotation system, and a chord of that kind can be absent from the graph. After the mesh is built, each boundary loop of at most `max_hole_edges` edges (default 100) is replaced by a face and ear-clipped. A longer boundary stays open.

A stored neighbor list is reused only when it is an exact search of the requested count. The kD-tree heap drops distance ties once it is full, so a prefix of a larger search is a different neighbor set. A rank distance, such as the one-ring radius or the `2/3` length cap, may be read from a longer search.

## Owl-little, one collapse, Euclidean, genus 0, reexpansion skipped

The coarse meshes are the same size. The times were measured while every reconstruction edge still sat in the contraction star: about 1.7 s on the graph path against 0.15 s for this cloud of 10k points. The queue now holds only the one-ring. That run has not been repeated. The mesh counts below are from the skipped-reexpansion comparison.

| | vertices | faces | boundaries | graph edges | face candidates |
| --- | --- | --- | --- | --- | --- |
| `hrsr_recon` | 5001 plus one triangle | 9754 | 11 | rebuilt nearest-neighbor graph | 101032 |
| `hrsr_recon_graph` | 5004 | 9754 | 10 | 105312 | 93840 |

## Reexpansion

Reexpansion is on by default. `skip_reexpansion=True` leaves the coarse mesh, which is the mesh in the table above. On the clouds tried so far that coarse mesh is the better surface. Reinsertion is what damages it.

Each coarse vertex sits at `v_bar`, the weighted average stored for that collapse. Reexpansion replays the log in reverse and finds the vertex by exact `double` equality with `v_bar`. The coarse normal is saved first. The survivor is moved to the active sample, the vertex is split, and the new vertex is placed at the latent sample. The original vertex stays at the active sample.

The link of both vertices is then pushed on a queue. For an interior edge the two outer edges of the quadrilateral are queued before a flip, and the flipped edge itself is not. The outside face is the one that contains neither split vertex. When both or neither do, it is the face whose third vertex is farther from `v_bar`. An interior edge flips only when that strictly raises the minimum of the six corner angles. Boundary edges are left in place. The walk stops after 200 flips.

The split gap minimizes a dihedral cost that is zero on a flat surface, plus a minimum-angle cost. The angle cost is the smallest corner's deviation from 60°, times the shortest edge. A corner under 3° adds `min_angle_threshold_penalty` (0.1) on top of that deviation. The score is still multiplied by the shortest edge, so a sliver whose short edge is tiny can remain cheap. `refinement_iterations` is unused by this repair.

The position map is updated only after a successful split, so a record that finds no split gap drops the rest of that vertex's collapse chain.
