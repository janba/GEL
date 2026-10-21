# Local rotation system reconstruction

LRsR builds a triangle mesh from a point cloud. The candidate triangles come from the radius-limited Delaunay graph. A tunnel-weighted spanning forest seeds the fill. Triangles are then added smallest score first, and each connected component is given one consistent orientation.

The C++ entry point is `HMesh::LRsR::reconstruct` in `src/GEL/HMesh/LRsR.h`. It clears `result`, writes the mesh, prints `LRsR iterations: N`, and returns `N`. An empty cloud, or one with more points than `int` can hold, prints `LRsR iterations: 0` and returns 0.

In Python the function is `hmesh.lrsr_recon(vertices, radius_factor=4.0)`. It returns only the `Manifold`. The iteration line is the C++ print. `src/demo/LRsR/lrsr.py` reconstructs an OBJ point cloud, or 400 points on the unit sphere, and draws the mesh with the same forest.

`LRsROpts::radius_factor` defaults to 4, the factor the 3DDelaunay demo uses for a point cloud.

## Candidates and the forest

The Delaunay radius is `radius_factor` times the median nearest-neighbor distance. A non-positive factor sets the radius to 0. `Geometry::delaunay_edges` supplies the candidate graph, one node per input point. The third vertex of a triangle on an edge is a Delaunay neighbor of both endpoints.

Normals are the unoriented adaptive estimate: `RSR::point_cloud_normals` with Euclidean distance and `NormalOrientation::None`. That is the same estimate as `hmesh.rsr_normals(..., use_Euclid_dist=True, orientation="none")`. The run prints `adaptive normal median k`. A missing or short normal array is replaced by zeros.

The edge normal of undirected edge `a--b`, with `a < b`, is the average of the two endpoint normals. When their dot product is negative, the higher-index endpoint is flipped before the average. An edge whose normal is missing, or whose length is at most `1e-12`, is not a seed.

The forest is `Geometry::minimum_spanning_forest` on `Geometry::tunnel_weight_graph`. An edge of Euclidean length `L` costs `L * (1 + |ê·n_a| + |ê·n_b|)`. An edge that runs along the normals costs more than one that stays in the tangent plane, so the forest prefers to stay on a sheet.

## Score and sides

Left of directed edge `a → b`, with `a < b`, is the side where `(b - a) × (c - a)` agrees with the edge normal. The cosine of that cross product with the normal must exceed `1e-8`. Right is the opposite sign. A corner this close to the plane of the edge and its normal is on neither side.

The score to minimize is perimeter times triangle area, divided by the area projected along the edge normal. The projected area is the shadow of the triangle on the plane perpendicular to that normal. A zero or non-finite projected area is not a candidate. The same three vertices can receive a different score from another edge, because the edge normal differs. Ties break by the smaller corner index, then the smaller edge endpoints, then the side.

While a triangle is being chosen, the left triangle is stored as `(a, b, corner)` and the right triangle as `(a, corner, b)`. Those two windings run in opposite directions along the edge, and each newly chosen triangle has the edge normal on its positive side. The orientation pass below may reverse a triangle afterwards. Neighboring edges are not forced to agree during the fill.

## The queue

One priority queue holds candidate triangles. Each entry is one corner proposed by one side of one edge, scored with that edge's projected-perimeter score. At most one entry per side is live. The side keeps the remaining corners in score order.

The queue starts with both open sides of every forest edge. A side with no finite positive-area corner pushes nothing. The smallest score is popped.

A popped triangle is rejected, and that side's next corner is pushed, when any of these holds:

- It would lie on a side of any of its three edges that already holds a triangle.
- It would be a third triangle on any of those edges.
- At any of its three vertices, its wedge overlaps a triangle already stored there.

The wedge test projects the incident edges onto the plane perpendicular to the vertex normal. A spoke lies in a wedge when it is a strictly positive combination of the other two projected edges, past `1e-8`. The wedge is at most 180 degrees. Flipping the vertex normal does not change the test. An angle that collapses in that plane does not block. A shared edge is the boundary and is allowed. The test runs both ways, so a new triangle that contains an existing one is rejected along with one that falls inside it.

Rejection only gets stricter: occupancy, face counts, and fans never loosen, so a rejected corner is not tried again on that side. If the proposing side is already occupied, or this exact triangle is already in the mesh, the side is closed and no worse corner is tried.

The two sides of one edge are separate queue entries. The better score is added first. The other side then sees that triangle, so overlapping wedges are rejected by the same fan test. An edge may still receive both sides when the wedges do not overlap and neither addition would be a third face.

Geometric filters that do not depend on the mesh — positive area, a cosine past `1e-8`, and a finite score — are applied once, when a side's corner list is built.

## Growth

Accepting a triangle marks the side it occupies on each of its three edges. Each of the other two edges that is not already under consideration is added to that set. The new triangle occupies one side of such an edge under that edge's own normal. Only the opposite side is queued. An edge already under consideration is left alone: a forest edge queued both sides at the start, and an edge opened earlier already queued its open side.

The queue runs until it is empty. There is no pass cap.

`N`, the value printed as `LRsR iterations: N` and returned by `reconstruct`, is how many accepted triangles were proposed by an edge outside the forest. A triangle proposed by a forest edge is not counted.

## Orientation

After the queue is empty, each connected component of the triangle mesh is oriented. Two triangles are neighbors when they share an edge. The earliest triangle in the component, the one with the smallest index in the face list, keeps its cyclic order. The walk is a depth-first flood fill on an explicit stack.

A neighbor is reversed when its cyclic order traverses the shared edge in the same direction as the triangle already accepted into the component. Reversal swaps the second and third vertices. The walk then continues from that neighbor. A triangle already visited is not reversed again. The test uses only the cyclic order of the three vertex indices. Normals are not consulted.

Consistent neighbors run in opposite directions along every shared edge. A component that cannot be oriented, because some later edge disagrees with a triangle already fixed, keeps the orientation assigned on the first visit.

## Mesh output

Each oriented triangle is added with `Manifold::add_face`. `stitch_mesh` then clusters boundary vertices within `1e-8`, and `cleanup` runs. Interior edges whose windings agree can stitch. Duplicate vertices can remain where a boundary does not meet within that radius.

The vertex-fan rules above keep overlapping wedges out of a point's link. They do not require that link to be a single chain. Several disjoint fans around one point are allowed.
