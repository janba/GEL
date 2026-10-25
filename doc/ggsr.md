# Gabriel graph surface reconstruction
_written by Grok based on an algorithm by Andreas_

GGSR builds a polygon mesh from a point cloud. The faces are cycles traced on the radius-limited Gabriel graph. Every pair of edges at a vertex seeds a walk. The walk follows the turn that finishes the shortest cycle, and each connected set of faces is given one consistent orientation.

The C++ entry point is `HMesh::GGSR::reconstruct` in `src/GEL/HMesh/GGSR.h`. It clears `result`, writes the mesh, prints `GGSR mst edges: M` and `GGSR faces: N`, and returns `N`. An empty cloud, or one with more points than `int` can hold, prints both lines as 0 and returns 0.

In Python the function is `hmesh.ggsr_recon(vertices, radius_factor=4.0, max_edges=10)`. It returns only the `Manifold`. The two count lines are the C++ print. `src/demo/GGSR/ggsr.py` reconstructs an OBJ point cloud, or 400 points on the unit sphere, and draws the mesh together with the forest. `--radius` sets the factor and `--edges` sets `max_edges`.

`GGSROpts::radius_factor` defaults to 4. `GGSROpts::max_edges` defaults to 10 and is the maximum number of edges on one traced face.

## Graph, normals, and the forest

The Gabriel radius is `radius_factor` times the median nearest-neighbor distance. A non-positive factor sets the radius to 0. With fewer than two points, or with no finite positive neighbor distance, the spacing used in that product is 1. An even count of nearest distances takes the average of the two central values. `Geometry::gabriel_edges` supplies the graph, one node per input point, in input order.

Normals are the unoriented adaptive estimate: `RSR::point_cloud_normals` with Euclidean distance and `NormalOrientation::None`. That is the same estimate as `hmesh.rsr_normals(..., use_Euclid_dist=True, orientation="none")`. The run prints `adaptive normal median k`. A missing or short normal array is replaced by zeros.

The forest is `Geometry::minimum_spanning_forest` on `Geometry::tunnel_weight_graph`. An edge of Euclidean length `L` costs `L * (1 + |ê·n_a| + |ê·n_b|)`. An edge that runs along the normals costs more than one that stays in the tangent plane. `M`, printed as `GGSR mst edges: M`, is the number of forest edges. Face finding does not use the forest. The demo draws it beside the mesh.

## Corners

A corner is an unordered pair of Gabriel edges that meet at a vertex. Every pair whose two directions are defined is a corner. A direction whose squared length is at most `1e-24`, or is not finite, is skipped. A vertex of degree less than 2 has no corner. Neighbors come out of the adjacency map in increasing id order, and each pair is stored in that order.

Both edges of a corner seed a walk. One walk leaves along the lower id. The other leaves along the higher id.

## The walk

The seed normal is the cross product of the departure edge with the other corner edge, normalized. The corner then lies to the left of travel. A cross product whose squared length is at most `1e-24`, or is not finite, skips that seed.

At the current vertex every other incident edge is a candidate. The previous edge is not. The inbound and outbound directions are projected onto the plane perpendicular to the running normal. A projection whose squared length is at most `1e-24`, or is not finite, is skipped. The score is the signed turn

```
atan2(n · (in × out), in · out)
```

with `in` and `out` already projected. A positive score turns left of travel. A score above `1e-9` counts as progress toward closing. A straight or right turn does not.

The candidates are ranked by three keys, in order. The smallest close length wins. An equal length keeps the greatest signed turn. An equal turn keeps the smaller neighbor id.

The close length is assigned as follows.

- The next vertex is already on this walk, the turn is positive, and the cycle from that earlier vertex through the current one has between 3 and `max_edges` edges. The close length is that cycle's edge count. When `kExclusiveCorners` is set, a taken corner at the earlier vertex removes this candidate.
- Otherwise the turn is positive, the next vertex is not on the walk, it is a Gabriel neighbor of the start, and taking it would leave a cycle of between 3 and `max_edges` edges. The close length is the number of edges already walked plus two. The face would then close on the following step.
- Otherwise the close length is `max_edges + 1`. The greatest signed turn decides.

A chosen next vertex that already lies on the walk ends it. The returned cycle is the part from that vertex through the current one. The vertices walked before that point are dropped, so the seed corner need not lie on the face. Returning to the start is this same rule: the start is the first vertex of the walk.

If no candidate survives, or the chosen corner at the current vertex is already claimed, the walk returns nothing. A repeated vertex that is not a positive close also ends the walk.

After a step that does not close, the running normal is replaced by the cross product of the two accepted edges at the current vertex. The cross product is flipped when its dot product with the old normal is negative, so the new normal stays in the same closed hemisphere. A cross product whose squared length is at most `1e-24`, or is not finite, leaves the normal unchanged.

The walk takes at most `max_edges` steps. When the graph has fewer nodes than `int` can hold, the limit is also capped by the number of nodes. `max_edges` below 3 yields no faces. A kept cycle is simple and has at least three edges.

## One corner, one cycle

`kExclusiveCorners` in `src/GEL/HMesh/GGSR.cpp` is currently `false`. While it is false, a corner may belong to several cycles. The checks and the claim are still in the tracer and return immediately.

When the flag is set, a finished cycle claims the corner at each of its vertices: the unordered pair of its two cycle neighbors there. A later walk does not start from a claimed corner. It stops when the corner it would step onto is already claimed, rather than taking another edge. The first completed walk keeps the corner.

## Keeping a cycle once

Each kept cycle is rotated so its least vertex id is first. If the second vertex is greater than the last, the rest of the cycle is reversed. The cycle then leaves its least vertex toward the smaller of its two cycle neighbors.

The cycles are sorted and identical sequences are dropped. `N`, printed as `GGSR faces: N` and returned by `reconstruct`, is the number of cycles left after that. It is counted from that list. Stitching does not change the printed number.

## Orientation

After duplicates are dropped, each connected set of faces is oriented. Two faces are neighbors when they share an edge. The earliest face in the set, the lexicographically smallest canonical cycle, keeps its order. The walk is a breadth-first flood fill.

A neighbor is reversed when its cyclic order traverses the shared edge in the same direction as the face already accepted. Reversal reverses the whole vertex list. The walk then continues from that neighbor. A face already visited is not reversed again. The test uses only the cyclic order of the vertex indices. Normals are not consulted.

Consistent neighbors run in opposite directions along a shared edge. A set that cannot be oriented, because some later edge disagrees with a face already fixed, keeps the orientation assigned on the first visit. An edge in three or more faces cannot give every pair opposite directions.

## Mesh output

Each oriented polygon is added with `Manifold::add_face`, using the input point positions. `stitch_mesh` then clusters boundary vertices within `1e-8`, and `cleanup` runs. Vertices that do not meet within that radius stay separate. An edge that would belong to a third face cannot be stitched into a manifold link, so a vertex can appear twice in a one-ring and `valid` can be false.

There is no larger-face discard and no cap on how many faces share an edge. Not every cycle of the Gabriel graph is found. A through-sheet edge, one that runs along the vertex normals, is an ordinary Gabriel edge. A corner that includes one still seeds a walk, and that walk can close a face through the sheet instead of along it.
