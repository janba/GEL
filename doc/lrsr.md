# Local rotation system reconstruction
_written by Grok based on an algorithm by Andreas_

LRsR builds a triangle mesh from a point cloud. The candidate triangles come from the radius-limited Delaunay graph. `kGrowthMode` in `src/GEL/HMesh/LRsR.cpp` selects which edges may propose triangles. It is `2`. Triangles are added smallest score first, and each connected component is given one consistent orientation.

The C++ entry point is `HMesh::LRsR::reconstruct` in `src/GEL/HMesh/LRsR.h`. It clears `result`, writes the mesh, prints `LRsR mode: M` and `LRsR iterations: N`, and returns `N`. An empty cloud, or one with more points than `int` can hold, prints mode `M` and `LRsR iterations: 0` and returns 0.

In Python the function is `hmesh.lrsr_recon(vertices, radius_factor=4.0)`. It returns only the `Manifold`. The mode and iteration lines are the C++ print. `src/demo/LRsR/lrsr.py` reconstructs an OBJ point cloud, or 400 points on the unit sphere, and draws the mesh together with the tunnel-weighted spanning forest of the same Delaunay graph. Modes 1 and 3 build that same forest and use its edges. Mode 2 does not.

`LRsROpts::radius_factor` defaults to 4, the factor the 3DDelaunay demo uses for a point cloud.

## Candidates

The Delaunay radius is `radius_factor` times the median nearest-neighbor distance. A non-positive factor sets the radius to 0. `Geometry::delaunay_edges` supplies the candidate graph, one node per input point. The third vertex of a triangle on an edge is a Delaunay neighbor of both endpoints.

Normals are the unoriented adaptive estimate: `RSR::point_cloud_normals` with Euclidean distance and `NormalOrientation::None`. That is the same estimate as `hmesh.rsr_normals(..., use_Euclid_dist=True, orientation="none")`. The run prints `adaptive normal median k`. A missing or short normal array is replaced by zeros.

The edge normal of undirected edge `a--b`, with `a < b`, is the average of the two endpoint normals. When their dot product is negative, the higher-index endpoint is flipped before the average. An edge whose normal is missing, or whose length is at most `1e-12`, does not propose a triangle. Both of its sides stay closed if the edge is opened by a triangle that uses it.

## Score and sides

Left of directed edge `a → b`, with `a < b`, is the side where `(b - a) × (c - a)` agrees with the edge normal. The cosine of that cross product with the normal must exceed `1e-8`. Right is the opposite sign. A corner this close to the plane of the edge and its normal is on neither side.

The score to minimize is perimeter times triangle area, divided by the area projected along the edge normal. The projected area is the shadow of the triangle on the plane perpendicular to that normal. A zero or non-finite projected area is not a candidate. The same three vertices can receive a different score from another edge, because the edge normal differs. Ties break by the smaller corner index, then the smaller edge endpoints, then the side.

While a triangle is being chosen, the left triangle is stored as `(a, b, corner)` and the right triangle as `(a, corner, b)`. Those two windings run in opposite directions along the edge, and each newly chosen triangle has the edge normal on its positive side. Mode 2 replaces that winding, once a parent face exists, with the order that traverses the shared edge the other way. Modes 1 and 3 keep the edge-normal winding until the orientation pass below.

## The queue

One priority queue holds candidate triangles. Each entry is one corner proposed by one side of one edge. Its score is that edge's projected-perimeter score times a factor for every face already across one of the triangle's three edges. At most one entry per side is live. The side keeps the remaining corners in score order.

The factor for one neighbor is `2 - cos(alpha)`. Alpha is the angle between the two face normals after the faces are given opposite directions along the shared edge, so the pair has one consistent orientation. Coplanar faces whose normals then agree use factor 1. Coplanar faces whose normals then oppose use factor 3. A fold lands between those. A triangle with no neighbor yet keeps factor 1. Several neighbors multiply.

The queue starts empty. What is pushed first depends on the mode, described below. Ties break by the smaller corner index, then the smaller edge endpoints, then the side. The smallest score is popped.

A neighbor can appear while an entry waits, and the factor then rises. Popping an entry recomputes the score. If the current score is higher, that entry is discarded and the same corner is pushed again at the current score. It is considered when that score reaches the front of the queue. An older copy of the same side is ignored.

A popped triangle is rejected, and that side's next corner is pushed, when any of these holds:

- It would lie on a side of any of its three edges that already holds a triangle.
- It would be a third triangle on any of those edges.
- At any of its three vertices, its wedge overlaps a triangle already stored there.

The wedge test projects the incident edges onto the plane perpendicular to the vertex normal. A spoke lies in a wedge when it is a strictly positive combination of the other two projected edges, past `1e-8`. The wedge is at most 180 degrees. Flipping the vertex normal does not change the test. An angle that collapses in that plane does not block. A shared edge is the boundary and is allowed. The test runs both ways, so a new triangle that contains an existing one is rejected along with one that falls inside it.

Rejection only gets stricter: occupancy, face counts, and fans never loosen, so a rejected corner is not tried again on that side. If the proposing side is already occupied, or this exact triangle is already in the mesh, the side is closed and no worse corner is tried.

The two sides of one edge are separate queue entries. The better score is added first. The other side then sees that triangle, so overlapping wedges are rejected by the same fan test. An edge may still receive both sides when the wedges do not overlap and neither addition would be a third face.

Geometric filters that do not depend on the mesh — positive area, a cosine past `1e-8`, and a finite score — are applied once, when a side's corner list is built.

## Growth

`kGrowthMode` is the switch.

Mode 1 keeps one queue. It starts with both sides of every edge of the tunnel-weighted spanning forest, built by `Geometry::minimum_spanning_forest` on `Geometry::tunnel_weight_graph`. An edge of Euclidean length `L` costs `L * (1 + |ê·n_a| + |ê·n_b|)`. Accepting a triangle also queues the free side of each of its other edges. A tree triangle and a triangle on a newly opened edge compete by score.

Mode 2 is the expanding front. The seed is the proposal with the lowest score: one corner, on one side of one Delaunay edge. The same three vertices can be scored by each of their three edges, and the winning proposal is the lowest of those scores. That triangle is emitted directly, with the edge-normal winding, and is not placed on the queue. Its three edges are the front. Each queues only the side the seed does not occupy. Later triangles queue the free side of an edge that is not already in the front. A sheet the seed does not reach stays empty. A side with no finite positive-area corner is not a candidate for the seed. If no such triangle exists, the mesh stays empty.

In mode 2 the winding is fixed when the triangle is queued. It traverses the proposing edge opposite the face already across that edge, the triangle that begot it. When the entry is popped, that winding is compared with every face already across any of its three edges. The same direction on any shared edge discards the triangle, and the next corner on that side is tried. A risen score is still pushed back before that test only when the winding still agrees. The check runs first: an inconsistent triangle is not put back.

Mode 3 uses the same forest as mode 1. Round one queues only the two sides of each tree edge. Accepting a triangle does not open its other edges, so a non-tree edge cannot propose until the tree queue is empty. Round two then queues the free side of each edge of the triangles already added, and from there continues like mode 2: each newly accepted triangle opens the free side of its edges that are not already in the front.

In every mode, accepting a triangle marks the side it occupies on each of its three edges. When an edge is opened, the new triangle occupies one side under that edge's own normal, and only the opposite side is queued. An edge already in the front is left alone. An edge with no usable normal, or one the new triangle does not sit clearly to the left or right of, closes both sides and queues nothing.

The queue runs until it is empty. There is no pass cap. Mode 3 runs that loop twice.

`N`, the value printed as `LRsR iterations: N` and returned by `reconstruct`, is how many triangles the queue accepted. Mode 2 does not count its first triangle. Modes 1 and 3 count every accepted triangle. The selected mode is printed as `LRsR mode: M`.

## Orientation

Mode 2 carries orientation with the front, as described above. The seed keeps its winding. Every later triangle is stored so that it already opposes its parent across the edge that queued it.

After the queue is empty, each connected component of the triangle mesh is oriented. Two triangles are neighbors when they share an edge. The earliest triangle in the component, the one with the smallest index in the face list, keeps its cyclic order. The walk is a depth-first flood fill on an explicit stack.

A neighbor is reversed when its cyclic order traverses the shared edge in the same direction as the triangle already accepted into the component. Reversal swaps the second and third vertices. The walk then continues from that neighbor. A triangle already visited is not reversed again. The test uses only the cyclic order of the three vertex indices. Normals are not consulted.

Consistent neighbors run in opposite directions along every shared edge. A component that cannot be oriented, because some later edge disagrees with a triangle already fixed, keeps the orientation assigned on the first visit.

## Mesh output

Each oriented triangle is added with `Manifold::add_face`. `stitch_mesh` then clusters boundary vertices within `1e-8`, and `cleanup` runs. Interior edges whose windings agree can stitch. Duplicate vertices can remain where a boundary does not meet within that radius.

The vertex-fan rules above keep overlapping wedges out of a point's link. They do not require that link to be a single chain. Several disjoint fans around one point are allowed.
