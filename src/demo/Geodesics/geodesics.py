"""Geodesics by the heat method of Crane, Weischedel, and Wardetzky.

The first click on the mesh is the source. The program integrates the heat
equation from that point, normalizes the heat gradient on each triangle, and
solves a Poisson equation for the divergence of that unit field. The solution
is the geodesic distance. Each later click traces the distance downhill, back
to the source, and adds that curve in the viewer.
"""

import math
from sys import argv

import colorsys
import numpy as np
import polyscope as ps
import polyscope.imgui as psim
from pygel3d import hmesh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

EPS = 1e-12
# Share of each step taken from the current triangle. The rest is the direction
# carried in, so a few degrees of gradient noise cannot bend a flat geodesic
# onto a mesh edge.
DIRECTION_MIX = 0.15


def barycentric(p0, p1, p2, point):
    """Barycentric coordinates of point in the plane of triangle p0, p1, p2."""
    v0 = p1 - p0
    v1 = p2 - p0
    v2 = point - p0
    d00 = np.dot(v0, v0)
    d01 = np.dot(v0, v1)
    d11 = np.dot(v1, v1)
    d20 = np.dot(v2, v0)
    d21 = np.dot(v2, v1)
    denom = d00 * d11 - d01 * d01
    if abs(denom) < 1e-20:
        return None
    b1 = (d11 * d20 - d01 * d21) / denom
    b2 = (d00 * d21 - d01 * d20) / denom
    return np.array([1.0 - b1 - b2, b1, b2])


def sanitize(bary):
    clipped = np.maximum(bary, 0.0)
    total = clipped.sum()
    if total <= 0.0:
        return np.full(3, 1.0 / 3.0)
    return clipped / total


def point_at(mesh, face, bary):
    verts = mesh["faces"][face]
    corners = mesh["positions"][verts]
    return bary[0] * corners[0] + bary[1] * corners[1] + bary[2] * corners[2]


def shortest_edge(mesh, face):
    verts = mesh["faces"][face]
    corners = mesh["positions"][verts]
    return float(min(
        np.linalg.norm(corners[1] - corners[0]),
        np.linalg.norm(corners[2] - corners[1]),
        np.linalg.norm(corners[0] - corners[2]),
    ))


def build_mesh(manifold):
    """Triangle corners and the face across each edge.

    Slot i of a face stores the neighbor across the edge opposite vertex i.
    """
    positions = np.array(manifold.positions(), dtype=np.float64, copy=True)
    faces = []
    compact = {}
    for face in manifold.faces():
        face = int(face)
        incident = [int(v) for v in manifold.circulate_face(face, "v")]
        if len(incident) != 3:
            continue
        compact[face] = len(faces)
        faces.append(incident)
    faces = np.asarray(faces, dtype=np.int64)
    neighbors = np.full((len(faces), 3), -1, dtype=np.int64)
    for face, index in compact.items():
        adjacent = [int(other) for other in manifold.circulate_face(face, "f")]
        if len(adjacent) != 3:
            continue
        for slot in range(3):
            other = adjacent[slot]
            if other in compact:
                neighbors[index, (slot + 1) % 3] = compact[other]
    vertex_faces = [[] for _ in range(len(positions))]
    for index, verts in enumerate(faces):
        for vertex in verts:
            vertex_faces[vertex].append(index)
    low = positions.min(axis=0)
    high = positions.max(axis=0)
    return {
        "positions": positions,
        "faces": faces,
        "neighbors": neighbors,
        "vertex_faces": vertex_faces,
        "diagonal": float(np.linalg.norm(high - low)),
    }


def hat_gradients(p0, p1, p2):
    """Gradients of the three hat functions, and the triangle area.

    The gradient of hat i is the inward perpendicular to the opposite edge,
    with length 1 / altitude.
    """
    normal = np.cross(p1 - p0, p2 - p0)
    twice_area = float(np.linalg.norm(normal))
    if twice_area < 1e-18:
        return None
    unit = normal / twice_area
    opposite = (p2 - p1, p0 - p2, p1 - p0)
    grads = tuple(np.cross(unit, edge) / twice_area for edge in opposite)
    return grads, 0.5 * twice_area


def build_operators(mesh):
    """Cotangent Laplacian, lumped mass, and the heat time step.

    The Laplacian is the stiffness matrix of the hat functions, so its action
    on a vertex scalar equals the integrated divergence of that scalar's
    gradient. The stored time step is the short time from the heat method,
    one mean edge length squared. heat_distance widens it when that kernel
    does not reach across the surface.
    """
    positions = mesh["positions"]
    faces = mesh["faces"]
    count = len(positions)
    rows = []
    cols = []
    values = []
    mass = np.zeros(count)
    hats = []
    lengths = []
    used = np.zeros(count, dtype=bool)
    for verts in faces:
        corners = positions[verts]
        lengths.extend((
            float(np.linalg.norm(corners[1] - corners[0])),
            float(np.linalg.norm(corners[2] - corners[1])),
            float(np.linalg.norm(corners[0] - corners[2])),
        ))
        local = hat_gradients(*corners)
        hats.append(local)
        if local is None:
            continue
        grads, area = local
        for vertex in verts:
            used[int(vertex)] = True
            mass[int(vertex)] += area / 3.0
        for a in range(3):
            for b in range(3):
                rows.append(int(verts[a]))
                cols.append(int(verts[b]))
                values.append(area * float(np.dot(grads[a], grads[b])))
    for vertex in range(count):
        if not used[vertex]:
            rows.append(vertex)
            cols.append(vertex)
            values.append(1.0)
            mass[vertex] = 1.0
    laplacian = coo_matrix((values, (rows, cols)), shape=(count, count)).tocsc()
    laplacian.sum_duplicates()
    mass[mass <= 0.0] = 1.0
    mean_edge = float(np.mean(lengths)) if lengths else 1.0
    mesh["hats"] = hats
    mesh["mass"] = mass
    mesh["mean_edge"] = mean_edge
    mesh["time"] = mean_edge * mean_edge
    return laplacian


def scalar_gradient(mesh, values, face):
    local = mesh["hats"][face]
    if local is None:
        return None
    grads, _ = local
    verts = mesh["faces"][face]
    return values[verts[0]] * grads[0] + values[verts[1]] * grads[1] + values[verts[2]] * grads[2]


def integrated_divergence(mesh, field):
    """Integrated divergence of one vector per face. Pairs with the Laplacian."""
    divergence = np.zeros(len(mesh["positions"]))
    for face, vector in enumerate(field):
        local = mesh["hats"][face]
        if local is None:
            continue
        grads, area = local
        verts = mesh["faces"][face]
        for slot in range(3):
            divergence[int(verts[slot])] += area * float(np.dot(grads[slot], vector))
    return divergence


def solve_poisson(laplacian, rhs, pin):
    """Solve L phi = rhs with phi[pin] = 0."""
    mask = np.ones(laplacian.shape[0], dtype=bool)
    mask[pin] = False
    keep = np.flatnonzero(mask)
    system = laplacian[keep][:, keep].tocsc()
    solution = spsolve(system, rhs[keep])
    if not np.all(np.isfinite(solution)):
        system = system + diags(np.full(len(keep), 1e-8))
        solution = spsolve(system, rhs[keep])
    phi = np.zeros(laplacian.shape[0])
    phi[keep] = solution
    return phi


def unit_field(mesh, heat):
    """Unit vectors pointing away from the heat peak, one per face."""
    field = []
    missing = 0
    for face in range(len(mesh["faces"])):
        gradient = scalar_gradient(mesh, heat, face)
        if gradient is None:
            field.append(np.zeros(3))
            missing += 1
            continue
        length = float(np.linalg.norm(gradient))
        if length < 1e-12:
            field.append(np.zeros(3))
            missing += 1
        else:
            field.append(-gradient / length)
    return field, missing


def shift_distance(mesh, distance, face, bary):
    """Move the distance so the source point is zero and values increase away from it."""
    distance = np.asarray(distance, dtype=np.float64).copy()
    distance -= float(np.dot(bary, distance[mesh["faces"][face]]))
    if float(np.mean(distance)) < 0.0:
        distance = -distance
    return distance


def distance_slope(mesh, distance):
    """Tenth percentile of |grad distance|. A finished field is near 1."""
    slopes = np.empty(len(mesh["faces"]))
    for face in range(len(mesh["faces"])):
        gradient = scalar_gradient(mesh, distance, face)
        slopes[face] = 0.0 if gradient is None else float(np.linalg.norm(gradient))
    if len(slopes) == 0:
        return 0.0
    return float(np.percentile(slopes, 10))


def heat_distance(mesh, laplacian, face, bary):
    """Geodesic distance from a point, by the three steps of the heat method.

    Integrate heat from a unit source at the point, replace the heat gradient
    by its unit direction away from the source, then solve a Poisson problem
    whose right-hand side is the divergence of that direction field. The time
    step starts at one squared mean edge length and doubles while the kernel
    still fails to reach the whole surface. The distance is kept once its
    gradient is near unit length on all but a few faces.
    """
    bary = sanitize(np.asarray(bary, dtype=np.float64))
    count = len(mesh["positions"])
    source = np.zeros(count)
    for weight, vertex in zip(bary, mesh["faces"][face]):
        source[int(vertex)] += float(weight)
    time = mesh["mean_edge"] ** 2
    if time <= 0.0:
        time = 1.0
    pin = int(mesh["faces"][face][int(np.argmax(bary))])
    distance = np.zeros(count)
    faces = max(1, len(mesh["faces"]))
    for _ in range(12):
        heat = spsolve((time * laplacian + diags(mesh["mass"])).tocsc(), source)
        field, missing = unit_field(mesh, heat)
        distance = solve_poisson(laplacian, integrated_divergence(mesh, field), pin)
        distance = shift_distance(mesh, distance, face, bary)
        if missing / faces <= 0.001 and distance_slope(mesh, distance) >= 0.7:
            break
        if time >= (mesh["mean_edge"] ** 2) * 4096.0:
            break
        time *= 2.0
    mesh["time"] = time
    return distance


def bary_rate(mesh, face, direction):
    """Change in barycentric coordinates per unit step along direction."""
    corners = mesh["positions"][mesh["faces"][face]]
    shifted = barycentric(*corners, corners[0] + direction)
    if shifted is None:
        return None
    return shifted - np.array([1.0, 0.0, 0.0])


def exit_along(mesh, face, bary, direction):
    """First boundary hit walking from bary along direction.

    Returns the opposite-vertex index of the exit edge, the step length, the
    exit barycentric coordinates, and the exit point. The step is zero when
    the point is already on an edge and direction leaves the face.
    """
    rate = bary_rate(mesh, face, direction)
    if rate is None:
        return -1, 0.0, None, None
    best_s = math.inf
    best_i = -1
    for slot in range(3):
        if rate[slot] >= -1e-14:
            continue
        step = -float(bary[slot]) / float(rate[slot])
        if step < best_s:
            best_s = step
            best_i = slot
    if best_i < 0 or not math.isfinite(best_s):
        return -1, 0.0, None, None
    best_s = max(0.0, best_s)
    exit_bary = sanitize(bary + best_s * rate)
    exit_bary[best_i] = 0.0
    exit_bary = sanitize(exit_bary)
    return best_i, best_s, exit_bary, point_at(mesh, face, exit_bary)


def sample(mesh, values, face, bary):
    return float(np.dot(bary, values[mesh["faces"][face]]))


def append_point(points, point):
    if np.linalg.norm(points[-1] - point) > 0.0:
        points.append(np.asarray(point, dtype=np.float64))


def clip_bary(bary):
    """Clamp to the triangle and renormalize. None when every weight is nonpositive."""
    clipped = np.maximum(np.asarray(bary, dtype=np.float64), 0.0)
    total = float(clipped.sum())
    if total <= 0.0:
        return None
    return clipped / total


def locate(bary):
    """Classify a point as a vertex, an edge, or a face interior.

    Returns the kind and the slot of the vertex, or of the edge's opposite
    vertex. The face-interior slot is unused.
    """
    slots = [slot for slot in range(3) if bary[slot] <= 1e-5]
    if len(slots) >= 2:
        return "vertex", int(np.argmax(bary))
    if len(slots) == 1:
        return "edge", slots[0]
    return "face", -1


def downhill(mesh, distance, face):
    """Unit direction of steepest decrease and its slope d(phi)/ds."""
    gradient = scalar_gradient(mesh, distance, face)
    if gradient is None:
        return None
    length = float(np.linalg.norm(gradient))
    if length < 1e-12:
        return None
    return -gradient / length, -length


def project_to_face(mesh, face, direction):
    """Unit projection of a direction onto the plane of a face."""
    corners = mesh["positions"][mesh["faces"][face]]
    normal = np.cross(corners[1] - corners[0], corners[2] - corners[0])
    length = float(np.linalg.norm(normal))
    if length < 1e-18:
        return None
    tangential = np.asarray(direction, dtype=np.float64) - normal * (
        float(np.dot(direction, normal)) / (length * length)
    )
    projected = float(np.linalg.norm(tangential))
    if projected < 1e-10:
        return None
    return tangential / projected


def mix_direction(carried, current):
    """Blend the incoming direction with the local downhill direction."""
    if current is None:
        return None if carried is None else np.asarray(carried, dtype=np.float64)
    if carried is None:
        return np.asarray(current, dtype=np.float64)
    blended = (1.0 - DIRECTION_MIX) * np.asarray(carried, dtype=np.float64)
    blended = blended + DIRECTION_MIX * np.asarray(current, dtype=np.float64)
    length = float(np.linalg.norm(blended))
    if length < 1e-12:
        return np.asarray(current, dtype=np.float64)
    return blended / length


def consensus_downhill(mesh, distance, vertex):
    """Mean of the unit downhill directions around a vertex."""
    total = np.zeros(3)
    for face in mesh["vertex_faces"][vertex]:
        descent = downhill(mesh, distance, face)
        if descent is not None:
            total += descent[0]
    length = float(np.linalg.norm(total))
    if length < 1e-10:
        return None
    return total / length


def sector_contains(rate, vertex_slot):
    """True when a direction at a vertex stays inside that face, edges included."""
    if rate is None or float(rate[vertex_slot]) >= -1e-8:
        return False
    for slot in range(3):
        if slot != vertex_slot and float(rate[slot]) < -1e-8:
            return False
    return True


def transfer_across(mesh, face, slot, bary):
    """Carry an edge point onto the neighboring face.

    The two edge weights are copied onto the same vertices. Nothing is
    reprojected, so a point in the middle of an edge cannot snap to a vertex.
    Returns the neighbor, the barycentric coordinates there, and the slot of
    the shared edge on the neighbor.
    """
    nxt = int(mesh["neighbors"][face, slot])
    if nxt < 0:
        return None
    verts = mesh["faces"][face]
    first = int(verts[(slot + 1) % 3])
    second = int(verts[(slot + 2) % 3])
    weights = {
        first: float(bary[(slot + 1) % 3]),
        second: float(bary[(slot + 2) % 3]),
    }
    carried = np.zeros(3)
    shared = -1
    matched = 0
    for index, vertex in enumerate(mesh["faces"][nxt]):
        vertex = int(vertex)
        if vertex in weights:
            carried[index] = weights[vertex]
            matched += 1
        else:
            shared = index
    if matched != 2 or shared < 0:
        return None
    carried = clip_bary(carried)
    if carried is None:
        return None
    return nxt, carried, shared


def accept_move(mesh, distance, face, bary, new_face, new_bary):
    """Keep a move only when it travels and the distance strictly falls."""
    new_bary = clip_bary(new_bary)
    if new_bary is None:
        return None
    moved = float(np.linalg.norm(
        point_at(mesh, new_face, new_bary) - point_at(mesh, face, bary)
    ))
    if moved <= 1e-8 * mesh["mean_edge"]:
        return None
    old_phi = sample(mesh, distance, face, bary)
    new_phi = sample(mesh, distance, new_face, new_bary)
    if new_phi >= old_phi - 1e-8 * mesh["mean_edge"]:
        return None
    return new_face, new_bary


def move_from_vertex(mesh, distance, face, slot):
    """Steepest feasible step away from a vertex.

    A face is feasible when its own downhill direction lies in that face.
    An edge is always feasible. The face wins a tie, so the trace enters a
    triangle whenever one contains the descent instead of snapping along a
    mesh edge. An edge is used when every incident downhill direction points
    out of its triangle, which is what used to end the curve early.
    """
    vertex = int(mesh["faces"][face][slot])
    phi = float(distance[vertex])
    candidates = []
    seen = set()
    for incident in mesh["vertex_faces"][vertex]:
        verts = mesh["faces"][incident]
        local = int(np.where(verts == vertex)[0][0])
        descent = downhill(mesh, distance, incident)
        if descent is not None:
            direction, slope = descent
            rate = bary_rate(mesh, incident, direction)
            if sector_contains(rate, local):
                start = np.zeros(3)
                start[local] = 1.0
                _edge, step, exit_bary, _point = exit_along(mesh, incident, start, direction)
                if exit_bary is not None and step > 1e-8 * mesh["mean_edge"]:
                    candidates.append((slope, 0, incident, exit_bary))
        for offset in (1, 2):
            other_slot = (local + offset) % 3
            other = int(verts[other_slot])
            key = (vertex, other) if vertex < other else (other, vertex)
            if key in seen:
                continue
            seen.add(key)
            length = float(np.linalg.norm(mesh["positions"][other] - mesh["positions"][vertex]))
            if length < 1e-15:
                continue
            slope = (float(distance[other]) - phi) / length
            if slope >= -1e-8:
                continue
            landed = np.zeros(3)
            landed[other_slot] = 1.0
            candidates.append((slope, 1, incident, landed))
    candidates.sort(key=lambda item: (item[0], item[1]))
    for _slope, _tie, new_face, new_bary in candidates:
        start = np.zeros(3)
        start[slot] = 1.0
        accepted = accept_move(mesh, distance, face, start, new_face, new_bary)
        if accepted is not None:
            return accepted
    return None


def move_from_edge(mesh, distance, face, bary, slot):
    """Cross into a face whose downhill direction enters it, or slide the edge.

    Sliding is the descent of a valley: both adjacent gradients point back
    onto the edge, so the piecewise-linear distance falls toward one endpoint.
    """
    options = [(face, bary, slot)]
    carried = transfer_across(mesh, face, slot, bary)
    if carried is not None:
        options.append(carried)
    candidates = []
    for other_face, other_bary, other_slot in options:
        descent = downhill(mesh, distance, other_face)
        if descent is None:
            continue
        direction, slope = descent
        rate = bary_rate(mesh, other_face, direction)
        if rate is None or float(rate[other_slot]) <= 1e-8:
            continue
        _edge, step, exit_bary, _point = exit_along(mesh, other_face, other_bary, direction)
        if exit_bary is None or step <= 1e-8 * mesh["mean_edge"]:
            continue
        candidates.append((slope, other_face, exit_bary))
    candidates.sort(key=lambda item: item[0])
    for _slope, new_face, new_bary in candidates:
        accepted = accept_move(mesh, distance, face, bary, new_face, new_bary)
        if accepted is not None:
            return accepted
    endpoints = ((slot + 1) % 3, (slot + 2) % 3)
    lower = min(endpoints, key=lambda index: float(distance[int(mesh["faces"][face][index])]))
    landed = np.zeros(3)
    landed[lower] = 1.0
    return accept_move(mesh, distance, face, bary, face, landed)


def move_inside(mesh, distance, face, bary):
    """Walk in a straight line to the first edge this face's downhill hits."""
    descent = downhill(mesh, distance, face)
    if descent is None:
        return None
    direction, _slope = descent
    _edge, step, exit_bary, _point = exit_along(mesh, face, bary, direction)
    if exit_bary is None or step <= 1e-8 * mesh["mean_edge"]:
        return None
    return accept_move(mesh, distance, face, bary, face, exit_bary)


def contains_point(mesh, face, point):
    """True when point lies in the triangle, boundary included."""
    corners = mesh["positions"][mesh["faces"][face]]
    bary = barycentric(*corners, point)
    if bary is None:
        return False
    return bool(np.min(bary) >= -1e-5)


def project_onto(mesh, face, point):
    verts = mesh["faces"][face]
    corners = mesh["positions"][verts]
    normal = np.cross(corners[1] - corners[0], corners[2] - corners[0])
    length = float(np.linalg.norm(normal))
    if length < 1e-18:
        return None
    unit = normal / length
    projected = point - unit * float(np.dot(point - corners[0], unit))
    bary = barycentric(*corners, projected)
    if bary is None:
        return None
    return sanitize(bary)


def direction_of_move(mesh, face, bary, new_face, new_bary, fallback):
    step = point_at(mesh, new_face, new_bary) - point_at(mesh, face, bary)
    length = float(np.linalg.norm(step))
    if length < 1e-15:
        return fallback
    return step / length


def descend_step(mesh, distance, face, bary, carried):
    """One downhill step, keeping most of the incoming direction.

    The carried direction is mixed with the local downhill direction and
    walked in a straight line to the next edge. Mixing stops the trace from
    snapping onto a mesh edge wherever the per-triangle gradient is only a
    few degrees off a straight geodesic. If that mixture does not descend,
    the step falls back to the steepest feasible move.
    """
    kind, slot = locate(bary)
    if kind == "vertex":
        vertex = int(mesh["faces"][face][slot])
        origin = np.zeros(3)
        origin[slot] = 1.0
        direction = mix_direction(carried, consensus_downhill(mesh, distance, vertex))
        best = None
        if direction is not None:
            for incident in mesh["vertex_faces"][vertex]:
                projected = project_to_face(mesh, incident, direction)
                if projected is None:
                    continue
                verts = mesh["faces"][incident]
                local = int(np.where(verts == vertex)[0][0])
                rate = bary_rate(mesh, incident, projected)
                if not sector_contains(rate, local):
                    continue
                start = np.zeros(3)
                start[local] = 1.0
                _edge, step, exit_bary, _point = exit_along(mesh, incident, start, projected)
                if exit_bary is None or step <= 1e-8 * mesh["mean_edge"]:
                    continue
                accepted = accept_move(mesh, distance, face, origin, incident, exit_bary)
                if accepted is None:
                    continue
                # Prefer the face in which this direction sits farthest from an edge.
                margin = min(float(rate[index]) for index in range(3) if index != local)
                if best is None or margin > best[0]:
                    best = (margin, accepted[0], accepted[1], projected)
        if best is not None:
            return best[1], best[2], best[3]
        moved = move_from_vertex(mesh, distance, face, slot)
        if moved is None:
            return None
        return moved[0], moved[1], direction_of_move(mesh, face, origin, moved[0], moved[1], direction)

    options = [(face, bary, -1 if kind == "face" else slot)]
    if kind == "edge":
        carried_edge = transfer_across(mesh, face, slot, bary)
        if carried_edge is not None:
            options.append(carried_edge)

    def consider(other_face, other_bary, other_slot, direction):
        if other_slot >= 0:
            rate = bary_rate(mesh, other_face, direction)
            if rate is None or float(rate[other_slot]) <= 1e-8:
                return
        _edge, step, exit_bary, _point = exit_along(mesh, other_face, other_bary, direction)
        if exit_bary is None or step <= 1e-8 * mesh["mean_edge"]:
            return
        accepted = accept_move(mesh, distance, face, bary, other_face, exit_bary)
        if accepted is None:
            return
        gradient = scalar_gradient(mesh, distance, other_face)
        slope = 0.0 if gradient is None else float(np.dot(gradient, direction))
        if best[0] is None or slope < best[0]:
            best[0] = slope
            best[1] = (accepted[0], accepted[1], direction)

    best = [None, None]
    for other_face, other_bary, other_slot in options:
        descent = downhill(mesh, distance, other_face)
        if descent is None:
            continue
        projected = None if carried is None else project_to_face(mesh, other_face, carried)
        mixed = mix_direction(projected, descent[0])
        if mixed is not None:
            consider(other_face, other_bary, other_slot, mixed)
        # A straight continuation still crosses a shallow valley when the
        # mixed direction has been bent back onto the edge.
        if projected is not None:
            gradient = scalar_gradient(mesh, distance, other_face)
            length = 0.0 if gradient is None else float(np.linalg.norm(gradient))
            if length > 1e-12 and float(np.dot(gradient, projected)) <= -0.5 * length:
                consider(other_face, other_bary, other_slot, projected)
    if best[1] is not None:
        return best[1]
    if kind == "edge":
        moved = move_from_edge(mesh, distance, face, bary, slot)
    else:
        moved = move_inside(mesh, distance, face, bary)
    if moved is None:
        return None
    return moved[0], moved[1], direction_of_move(mesh, face, bary, moved[0], moved[1], carried)


def closest_on_triangle(p0, p1, p2, point):
    """Closest point of a triangle to a query, including the boundary."""
    ab = p1 - p0
    ac = p2 - p0
    ap = point - p0
    d1 = float(np.dot(ab, ap))
    d2 = float(np.dot(ac, ap))
    if d1 <= 0.0 and d2 <= 0.0:
        return np.array(p0, dtype=np.float64, copy=True)
    bp = point - p1
    d3 = float(np.dot(ab, bp))
    d4 = float(np.dot(ac, bp))
    if d3 >= 0.0 and d4 <= d3:
        return np.array(p1, dtype=np.float64, copy=True)
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        span = d1 - d3
        if span <= 1e-18:
            return np.array(p0, dtype=np.float64, copy=True)
        v = d1 / span
        return p0 + v * ab
    cp = point - p2
    d5 = float(np.dot(ab, cp))
    d6 = float(np.dot(ac, cp))
    if d6 >= 0.0 and d5 <= d6:
        return np.array(p2, dtype=np.float64, copy=True)
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        span = d2 - d6
        if span <= 1e-18:
            return np.array(p0, dtype=np.float64, copy=True)
        w = d2 / span
        return p0 + w * ac
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        span = (d4 - d3) + (d5 - d6)
        if span <= 1e-18:
            return np.array(p1, dtype=np.float64, copy=True)
        w = (d4 - d3) / span
        return p1 + w * (p2 - p1)
    denom = va + vb + vc
    if denom <= 1e-18:
        choices = (p0, p1, p2)
        return np.array(min(choices, key=lambda vertex: float(np.dot(vertex - point, vertex - point))), copy=True)
    v = vb / denom
    w = vc / denom
    return p0 + ab * v + ac * w


def project_near(mesh, hint, point, rings=3, until_hit=False):
    """Closest point on the mesh in a neighborhood of hint.

    The result is a point of a visited triangle. It is never the raw query
    when that query sits off the surface. With until_hit, the search stops at
    the first face the point already lies on, so a sample several triangles
    away from the hint is not clamped to the wrong face.
    """
    point = np.asarray(point, dtype=np.float64)
    hint = int(hint)
    corners = mesh["positions"][mesh["faces"][hint]]
    best_point = closest_on_triangle(corners[0], corners[1], corners[2], point)
    best_dist = float(np.dot(best_point - point, best_point - point))
    best_face = hint
    tolerance = (1e-6 * max(mesh["diagonal"], 1.0)) ** 2
    if until_hit and best_dist <= tolerance:
        return best_face, best_point
    pending = []
    seen = {hint}
    for slot in range(3):
        other = int(mesh["neighbors"][hint, slot])
        if other >= 0 and other not in seen:
            seen.add(other)
            pending.append(other)
    for _ in range(max(0, int(rings))):
        if not pending:
            break
        for face in pending:
            corners = mesh["positions"][mesh["faces"][face]]
            on_face = closest_on_triangle(corners[0], corners[1], corners[2], point)
            delta = on_face - point
            dist = float(np.dot(delta, delta))
            if dist < best_dist:
                best_dist = dist
                best_face = face
                best_point = on_face
        if until_hit and best_dist <= tolerance:
            break
        nxt = []
        for face in pending:
            for slot in range(3):
                other = int(mesh["neighbors"][face, slot])
                if other >= 0 and other not in seen:
                    seen.add(other)
                    nxt.append(other)
        pending = nxt
    return best_face, best_point


def face_normal(mesh, face):
    """Unit normal of a triangle, or None when the triangle is degenerate."""
    corners = mesh["positions"][mesh["faces"][int(face)]]
    normal = np.cross(corners[1] - corners[0], corners[2] - corners[0])
    length = float(np.linalg.norm(normal))
    if length < 1e-18:
        return None
    return normal / length


def reproject(mesh, hint, point, rings=4, until_hit=False):
    """Snap point onto a triangle near hint.

    project_near already returns a point of a visited triangle. The snap is
    repeated here so a smoothing step cannot keep a point that still lies off
    that triangle.
    """
    face, on_face = project_near(mesh, hint, point, rings=rings, until_hit=until_hit)
    corners = mesh["positions"][mesh["faces"][int(face)]]
    snapped = closest_on_triangle(corners[0], corners[1], corners[2], np.asarray(on_face, dtype=np.float64))
    gap = float(np.linalg.norm(snapped - on_face))
    if gap > 1e-8 * max(mesh["mean_edge"], 1e-6):
        snapped = closest_on_triangle(corners[0], corners[1], corners[2], np.asarray(point, dtype=np.float64))
    return int(face), np.asarray(snapped, dtype=np.float64)


def tangential_reproject(mesh, face, point, target):
    """Move point toward target along the face, then snap back onto the mesh.

    The neighbor average of a curve on a curved surface points into the
    volume. Only the part lying in the tangent plane is kept, and that point
    is reprojected so the sample stays on a triangle.
    """
    point = np.asarray(point, dtype=np.float64)
    offset = np.asarray(target, dtype=np.float64) - point
    normal = face_normal(mesh, face)
    if normal is None:
        return reproject(mesh, face, point + offset, rings=3)
    slide = offset - normal * float(np.dot(offset, normal))
    proposed = point + slide
    reach = float(np.linalg.norm(slide))
    rings = max(3, int(math.ceil(reach / mesh["mean_edge"])) + 2)
    return reproject(mesh, face, proposed, rings=min(rings, 8))


def place_on_mesh(mesh, points, hint):
    """Snap each sample to the mesh, walking from the previous face."""
    face = int(hint)
    placed = []
    faces = []
    for point in points:
        face, on_face = reproject(mesh, face, point, rings=48, until_hit=True)
        placed.append(on_face)
        faces.append(face)
    placed[0] = np.asarray(points[0], dtype=np.float64)
    placed[-1] = np.asarray(points[-1], dtype=np.float64)
    return placed, faces


def face_holding(mesh, hint, point):
    """A face that contains point, searching outward from hint."""
    hint = int(hint)
    if contains_point(mesh, hint, point):
        return hint
    face, _on_face = project_near(mesh, hint, point, rings=8, until_hit=True)
    return int(face)


def _distance_to_chord(point, start, end):
    """Distance from point to the segment start-end."""
    direction = end - start
    length2 = float(np.dot(direction, direction))
    if length2 <= 1e-24:
        return float(np.linalg.norm(point - start))
    along = float(np.dot(point - start, direction) / length2)
    along = min(max(along, 0.0), 1.0)
    foot = start + along * direction
    return float(np.linalg.norm(point - foot))


def _sample_along(mesh, points, face_ids, cumulative, target, carried, step):
    """Point of the polyline at an arc length, snapped onto the mesh."""
    segment = 0
    while segment < len(points) - 2 and cumulative[segment + 1] < target - 1e-15:
        segment += 1
    span = cumulative[segment + 1] - cumulative[segment]
    alpha = 0.0 if span <= 1e-18 else (target - cumulative[segment]) / span
    alpha = min(max(alpha, 0.0), 1.0)
    guess = (1.0 - alpha) * points[segment] + alpha * points[segment + 1]
    # The sample can sit halfway along a much longer input segment, so the
    # search has to reach that far from the nearer end.
    if alpha <= 0.5:
        segment_hint = int(face_ids[segment])
        hop = alpha * span
    else:
        segment_hint = int(face_ids[segment + 1])
        hop = (1.0 - alpha) * span
    rings = min(48, max(4, int(math.ceil(max(hop, step) / mesh["mean_edge"])) + 4))
    carried_face, carried_point = reproject(mesh, carried, guess, rings=rings)
    segment_face, segment_point = reproject(mesh, segment_hint, guess, rings=rings)
    if float(np.dot(segment_point - guess, segment_point - guess)) < float(
        np.dot(carried_point - guess, carried_point - guess)
    ):
        return segment_face, segment_point
    return carried_face, carried_point


def _bend_to_keep(points, face_ids, cumulative, origin, target, start, end, gap_limit):
    """An input sample the chord start-end would cut away from the polyline.

    Dropping every sample between two stations is what turns a bend into a
    straight shortcut. Once that shortcut exists it is shorter than the path
    on the surface, so the next pass pulls harder and deletes still more
    samples. The sample farthest from the chord is kept when that distance
    exceeds the gap allowed of a chord.
    """
    best = None
    best_gap = gap_limit
    for index in range(1, len(points) - 1):
        arc = cumulative[index]
        if arc <= origin + 1e-9 or arc >= target - 1e-9:
            continue
        gap = _distance_to_chord(points[index], start, end)
        if gap > best_gap:
            best_gap = gap
            best = index
    if best is None:
        return None
    return (
        np.array(points[best], dtype=np.float64, copy=True),
        int(face_ids[best]),
        cumulative[best],
    )


def resample_curve(mesh, points, face_ids, spacing):
    """Place samples along the polyline, keeping a bend that a chord would cut.

    Averaging neighbors stops changing the length while clusters of samples
    are still packed together, so each shortening pass starts by evening the
    spacing. A sample is kept when the chord that would replace it passes more
    than a thousandth of the mean edge away from that sample. Those samples are
    what hold the curve on the surface; removing them leaves a straight
    shortcut through the mesh, and the shortcut is shorter, so the pull gets
    worse. Endpoints stay fixed. The spacing is widened when a curve is so
    long that equal steps would be finer than a four-hundredth of its length.
    """
    points = [np.asarray(point, dtype=np.float64) for point in points]
    if len(points) < 2:
        return points, [int(face) for face in face_ids]
    cumulative = [0.0]
    for start, end in zip(points, points[1:]):
        cumulative.append(cumulative[-1] + float(np.linalg.norm(end - start)))
    total = cumulative[-1]
    if total <= 1e-12:
        return (
            [np.array(points[0], copy=True), np.array(points[-1], copy=True)],
            [int(face_ids[0]), int(face_ids[-1])],
        )
    spacing = max(float(spacing), total / 400.0)
    count = max(1, int(round(total / spacing)))
    if len(points) >= 3:
        count = max(count, 2)
    step = total / count
    gap_limit = 0.001 * max(float(mesh["mean_edge"]), 1e-8)
    new_points = [np.array(points[0], dtype=np.float64, copy=True)]
    new_faces = [int(face_ids[0])]
    origin = 0.0
    carried = int(face_ids[0])
    # Each pass emits one kept bend sample or one evenly spaced sample.
    for _ in range(len(points) + count + 2):
        if origin >= total - 1e-9:
            break
        target = min(total, origin + step)
        face, on_face = _sample_along(mesh, points, face_ids, cumulative, target, carried, step)
        held = _bend_to_keep(
            points, face_ids, cumulative, origin, target, new_points[-1], on_face, gap_limit
        )
        if held is not None and held[2] > origin + 1e-9 and float(np.linalg.norm(held[0] - new_points[-1])) > 1e-12:
            new_points.append(held[0])
            new_faces.append(held[1])
            carried = held[1]
            origin = held[2]
            continue
        new_points.append(np.asarray(on_face, dtype=np.float64))
        new_faces.append(int(face))
        carried = int(face)
        origin = target
    new_points[0] = np.array(points[0], dtype=np.float64, copy=True)
    new_faces[0] = face_holding(mesh, face_ids[0], new_points[0])
    if float(np.linalg.norm(new_points[-1] - points[-1])) > 1e-12:
        new_points.append(np.array(points[-1], dtype=np.float64, copy=True))
        new_faces.append(int(face_ids[-1]))
    else:
        new_points[-1] = np.array(points[-1], dtype=np.float64, copy=True)
    new_faces[-1] = face_holding(mesh, face_ids[-1], new_points[-1])
    return new_points, new_faces


def _lay_triangle(p0, p1, p2):
    """Place a triangle in the plane with its real edge lengths."""
    edge = p1 - p0
    length = float(np.linalg.norm(edge))
    if length < 1e-18:
        return None
    to_apex = p2 - p0
    along = float(np.dot(to_apex, edge) / length)
    height = math.sqrt(max(float(np.dot(to_apex, to_apex)) - along * along, 0.0))
    return np.array([[0.0, 0.0], [length, 0.0], [along, height]], dtype=np.float64)


def _apex_across(edge_start, edge_end, dist_start, dist_end, other_side):
    """Third vertex of a triangle laid on the opposite side of an edge."""
    edge = edge_end - edge_start
    length = float(np.linalg.norm(edge))
    if length < 1e-18:
        return None
    direction = edge / length
    along = (dist_start * dist_start - dist_end * dist_end + length * length) / (2.0 * length)
    height = math.sqrt(max(dist_start * dist_start - along * along, 0.0))
    normal = np.array([-direction[1], direction[0]])
    side = float(direction[0] * (other_side[1] - edge_start[1]) - direction[1] * (other_side[0] - edge_start[0]))
    if side > 0.0:
        normal = -normal
    return edge_start + along * direction + height * normal


def _unfold(mesh, seeds, rings, reach, anchors):
    """Isometric layout of a neighborhood. Faces outside `reach` of every anchor are left out."""
    positions = mesh["positions"]
    faces = mesh["faces"]
    embedding = {}
    pending = []
    seen = set()
    for seed in seeds:
        seed = int(seed)
        laid = _lay_triangle(*positions[faces[seed]])
        if laid is None:
            continue
        embedding[seed] = laid
        pending.append(seed)
        seen.add(seed)
        break
    rings = max(1, int(rings))
    for _ in range(rings + 1):
        if not pending:
            break
        nxt = []
        for face in list(pending):
            if face not in embedding:
                continue
            emb = embedding[face]
            for slot in range(3):
                other = int(mesh["neighbors"][face, slot])
                if other < 0 or other in seen:
                    continue
                centroid = np.mean(positions[faces[other]], axis=0)
                if all(float(np.linalg.norm(centroid - anchor)) > reach for anchor in anchors):
                    continue
                start_slot = (slot + 1) % 3
                end_slot = (slot + 2) % 3
                other_verts = faces[other]
                shared = {int(faces[face, start_slot]), int(faces[face, end_slot])}
                new_slot = next((slot_index for slot_index in range(3) if int(other_verts[slot_index]) not in shared), None)
                if new_slot is None:
                    continue
                new_vertex = int(other_verts[new_slot])
                apex = _apex_across(
                    emb[start_slot],
                    emb[end_slot],
                    float(np.linalg.norm(positions[new_vertex] - positions[faces[face, start_slot]])),
                    float(np.linalg.norm(positions[new_vertex] - positions[faces[face, end_slot]])),
                    emb[slot],
                )
                if apex is None:
                    continue
                other_emb = np.zeros((3, 2))
                for slot_index in range(3):
                    vertex = int(other_verts[slot_index])
                    if vertex == int(faces[face, start_slot]):
                        other_emb[slot_index] = emb[start_slot]
                    elif vertex == int(faces[face, end_slot]):
                        other_emb[slot_index] = emb[end_slot]
                    else:
                        other_emb[slot_index] = apex
                embedding[other] = other_emb
                seen.add(other)
                nxt.append(other)
        pending = nxt
    return embedding


def _point_in_embedding(embedding, point):
    """Face and barycentric coordinates of a planar point, if it lies in the layout."""
    query = np.array([point[0], point[1], 0.0])
    best = None
    best_margin = -1e9
    for face, tri in embedding.items():
        bary = barycentric(
            np.array([tri[0, 0], tri[0, 1], 0.0]),
            np.array([tri[1, 0], tri[1, 1], 0.0]),
            np.array([tri[2, 0], tri[2, 1], 0.0]),
            query,
        )
        if bary is None:
            continue
        margin = float(np.min(bary))
        if margin > best_margin:
            best_margin = margin
            best = (int(face), bary)
    if best is None or best_margin < -1e-5:
        return None
    return best[0], sanitize(best[1])


def _bary_on_face(mesh, face, point):
    corners = mesh["positions"][mesh["faces"][face]]
    bary = barycentric(*corners, np.asarray(point, dtype=np.float64))
    if bary is None or float(np.min(bary)) < -1e-4:
        return None
    return sanitize(bary)


def _to_plane(tri, bary):
    return bary[0] * tri[0] + bary[1] * tri[1] + bary[2] * tri[2]


def geodesic_midpoint(mesh, face_a, point_a, face_c, point_c, hint_face, hint_point):
    """Midpoint of the short geodesic from point_a to point_c near hint_point.

    The neighborhood of the three samples is unfolded into the plane. The
    midpoint of the straight line is used when that line stays inside the
    unfolding. Otherwise there is no midpoint on the surface to follow. The
    Euclidean midpoint lies in the volume, and moving a sample toward it pulls
    the curve inside the mesh.
    """
    point_a = np.asarray(point_a, dtype=np.float64)
    point_c = np.asarray(point_c, dtype=np.float64)
    hint_point = np.asarray(hint_point, dtype=np.float64)
    chord = float(np.linalg.norm(point_c - point_a))
    if chord < 1e-12:
        return None
    face_a = int(face_a)
    face_c = int(face_c)
    hint_face = int(hint_face)
    if face_a == face_c:
        bary_a = _bary_on_face(mesh, face_a, point_a)
        bary_c = _bary_on_face(mesh, face_c, point_c)
        if bary_a is not None and bary_c is not None:
            mid_bary = sanitize(0.5 * (bary_a + bary_c))
            return face_a, point_at(mesh, face_a, mid_bary)
    reach = max(
        chord,
        float(np.linalg.norm(hint_point - point_a)),
        float(np.linalg.norm(hint_point - point_c)),
        mesh["mean_edge"],
    ) + 2.0 * mesh["mean_edge"]
    rings = int(math.ceil(reach / mesh["mean_edge"])) + 1
    rings = max(2, min(rings, 8))
    embedding = _unfold(mesh, (face_a, hint_face, face_c), rings, reach, (point_a, point_c, hint_point))
    bary_a = _bary_on_face(mesh, face_a, point_a) if face_a in embedding else None
    bary_c = _bary_on_face(mesh, face_c, point_c) if face_c in embedding else None
    if bary_a is not None and bary_c is not None:
        start = _to_plane(embedding[face_a], bary_a)
        end = _to_plane(embedding[face_c], bary_c)
        if all(_point_in_embedding(embedding, (1.0 - t) * start + t * end) is not None for t in (0.25, 0.5, 0.75)):
            located = _point_in_embedding(embedding, 0.5 * (start + end))
            if located is not None:
                face, bary = located
                return face, point_at(mesh, face, bary)
    # The straight line left the unfolding. An interior midpoint would shorten
    # the polyline by cutting through the mesh, and the straight chord would
    # then lose the samples that were holding the curve out.
    return None


def straighten_sweep(mesh, points, face_ids):
    """Move each interior sample toward the geodesic between its neighbors.

    The target has to lie on the surface. The step keeps only the tangential
    part of that displacement and reprojects onto the mesh. The sweep uses the
    updated position of the previous sample. A move that does not shorten the
    two adjacent segments is left out, so the polyline never grows.
    """
    points = [np.array(point, dtype=np.float64, copy=True) for point in points]
    faces = [int(face) for face in face_ids]
    tolerance = 1e-8 * mesh["mean_edge"]
    for index in range(1, len(points) - 1):
        located = geodesic_midpoint(
            mesh,
            faces[index - 1], points[index - 1],
            faces[index + 1], points[index + 1],
            faces[index], points[index],
        )
        if located is None:
            continue
        _face, target = located
        moved = tangential_reproject(mesh, faces[index], points[index], target)
        if moved is None:
            continue
        face, point = moved
        span = max(
            float(np.linalg.norm(points[index] - points[index - 1])),
            float(np.linalg.norm(points[index + 1] - points[index])),
        )
        if float(np.linalg.norm(point - points[index])) > span + mesh["mean_edge"]:
            continue
        before = (
            float(np.linalg.norm(points[index] - points[index - 1]))
            + float(np.linalg.norm(points[index + 1] - points[index]))
        )
        after = (
            float(np.linalg.norm(point - points[index - 1]))
            + float(np.linalg.norm(points[index + 1] - point))
        )
        if after > before + tolerance:
            continue
        points[index] = np.asarray(point, dtype=np.float64)
        faces[index] = int(face)
    return points, faces


def _best_chord_split(mesh, face_a, point_a, face_b, point_b):
    """Where a straight segment leaves the surface path, and the sample to insert.

    The short geodesic between the endpoints is unfolded. Samples of that
    geodesic that sit farthest from the chord are the bend the straight
    segment cuts through. The sample is kept only when it projects into the
    open middle of the chord, so the split cannot collapse onto an endpoint
    or jump to another sheet. When the unfolding does not cover the segment,
    the closest point of the mesh to the chord is used instead.
    """
    point_a = np.asarray(point_a, dtype=np.float64)
    point_b = np.asarray(point_b, dtype=np.float64)
    direction = point_b - point_a
    length2 = float(np.dot(direction, direction))
    if length2 <= 1e-24:
        return None
    length = math.sqrt(length2)
    best = None

    def consider(face, point):
        nonlocal best
        point = np.asarray(point, dtype=np.float64)
        along = float(np.dot(point - point_a, direction) / length2)
        if along < 0.05 or along > 0.95:
            return
        foot = point_a + along * direction
        gap = float(np.linalg.norm(point - foot))
        if gap > 0.75 * length:
            return
        if best is None or gap > best[2]:
            best = (int(face), point, gap)

    face_a = int(face_a)
    face_b = int(face_b)
    covered = False
    if face_a == face_b:
        bary_a = _bary_on_face(mesh, face_a, point_a)
        bary_b = _bary_on_face(mesh, face_b, point_b)
        if bary_a is not None and bary_b is not None:
            covered = True
            for station in (0.25, 0.5, 0.75):
                bary = sanitize((1.0 - station) * bary_a + station * bary_b)
                consider(face_a, point_at(mesh, face_a, bary))
    else:
        reach = max(length, mesh["mean_edge"]) + 2.0 * mesh["mean_edge"]
        rings = max(2, min(int(math.ceil(reach / mesh["mean_edge"])) + 1, 8))
        embedding = _unfold(mesh, (face_a, face_b), rings, reach, (point_a, point_b))
        bary_a = _bary_on_face(mesh, face_a, point_a) if face_a in embedding else None
        bary_b = _bary_on_face(mesh, face_b, point_b) if face_b in embedding else None
        if bary_a is not None and bary_b is not None:
            start = _to_plane(embedding[face_a], bary_a)
            end = _to_plane(embedding[face_b], bary_b)
            stations = (0.0625, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.9375)
            located = [_point_in_embedding(embedding, (1.0 - station) * start + station * end) for station in stations]
            if all(item is not None for item in located):
                covered = True
                for item in located:
                    face, bary = item
                    consider(face, point_at(mesh, face, bary))
    if covered:
        return best
    for station in (0.5, 0.25, 0.75, 0.375, 0.625):
        guess = (1.0 - station) * point_a + station * point_b
        hint = face_a if station <= 0.5 else face_b
        hop = min(station, 1.0 - station) * length
        rings = min(8, max(3, int(math.ceil(hop / mesh["mean_edge"])) + 3))
        face, on_mesh = reproject(mesh, hint, guess, rings=rings)
        consider(face, on_mesh)
    return best


def refine_high_curvature(mesh, points, face_ids):
    """Split chords that cut through a bend.

    Every sample already lies on a triangle. The drawn segment between two
    samples is straight, so it passes through the interior or through empty
    space wherever the surface bends. A segment gains an on-surface sample at
    the place the chord leaves the mesh the most, and the pieces are checked
    again. A chord closer than a thousandth of the mean edge is left as it is,
    so a flat stretch keeps its spacing. Splits stop at a five-hundredth of
    the mean edge. The endpoints stay fixed.
    """
    gap_limit = 0.001 * max(float(mesh["mean_edge"]), 1e-8)
    min_length = 0.002 * max(float(mesh["mean_edge"]), 1e-8)
    points = [np.array(point, dtype=np.float64, copy=True) for point in points]
    faces = [int(face) for face in face_ids]
    if len(points) != len(faces) or len(points) < 2:
        return points
    span = max(curve_length(points), float(mesh["mean_edge"]))
    cap = min(8000, max(len(points) + 4, int(math.ceil(span / min_length)) + 2))
    index = 0
    while index < len(points) - 1 and len(points) < cap:
        start = points[index]
        end = points[index + 1]
        length = float(np.linalg.norm(end - start))
        if length < 2.0 * min_length:
            index += 1
            continue
        split = _best_chord_split(mesh, faces[index], start, faces[index + 1], end)
        if split is None or split[2] <= gap_limit or split[2] > 0.75 * length:
            index += 1
            continue
        face, point, _gap = split
        if (
            float(np.linalg.norm(point - start)) < min_length
            or float(np.linalg.norm(end - point)) < min_length
        ):
            index += 1
            continue
        points.insert(index + 1, point)
        faces.insert(index + 1, face)
    return points


def shorten_curve(mesh, points, hint):
    """Shorten a polyline and even out its sample spacing.

    Each pass resamples the curve at equal steps, then moves every interior
    sample toward the geodesic between its neighbors. The move stays in the
    tangent plane and is reprojected onto a triangle, so the curve is not
    pulled into the volume. Coarse steps run first, so a bow that makes two
    geodesics cross can be cut off, and finer steps then take out the kinks.
    The endpoints stay fixed. At most ten passes run at a given spacing. A
    spacing stops early when a pass lengthens the curve, when the length drops
    by less than a ten-thousandth of the curve, or when a pass removes less
    than 15 percent of what the previous pass removed. Resampling keeps a
    sample when dropping it would leave a chord through the mesh, so a shortcut
    cannot shed the points that hold the curve out and then grow by being
    shorter. The result is resampled once more at half the mean edge length.
    Chords that still cut through a bend are then split, so high curvature is
    sampled finely enough that a straight segment stays on the surface.
    """
    if len(points) < 2:
        return [np.asarray(point, dtype=np.float64) for point in points]
    points, face_ids = place_on_mesh(mesh, points, hint)
    if len(points) >= 3:
        span = max(curve_length(points), mesh["mean_edge"])
        passes = 0
        for multiple, cap in ((4.0, 3), (2.0, 3), (1.0, 10)):
            if passes >= 10:
                break
            spacing = mesh["mean_edge"] * multiple
            previous_drop = None
            used = 0
            while passes < 10 and used < cap:
                length = curve_length(points)
                resampled, resampled_faces = resample_curve(mesh, points, face_ids, spacing)
                resampled_length = curve_length(resampled)
                if resampled_length <= length + 1e-8 * span:
                    points, face_ids = resampled, resampled_faces
                    length = resampled_length
                points, face_ids = straighten_sweep(mesh, points, face_ids)
                new_length = curve_length(points)
                drop = length - new_length
                passes += 1
                used += 1
                if drop <= 1e-4 * span:
                    break
                if previous_drop is not None and drop < 0.15 * previous_drop:
                    break
                previous_drop = drop
        points, face_ids = resample_curve(mesh, points, face_ids, 0.5 * mesh["mean_edge"])
    return refine_high_curvature(mesh, points, face_ids)


def trace_geodesic(mesh, distance, source_face, source_bary, start_face, start_bary):
    """Follow the distance downhill from a click to the source.

    The distance is linear on each triangle. The trace walks in a straight
    line across each one, carrying its direction forward and blending in a
    little of the local downhill direction. On an edge that neither side
    descends into, it slides toward the lower endpoint. At a vertex where no
    face contains the carried direction, it takes the steepest descending
    edge, so the curve does not die while a way downhill remains. The path
    joins the source once the source lies in the current triangle. The
    resulting polyline is then resampled at even steps and straightened. Each
    smoothing step stays in the tangent plane and is reprojected onto the mesh,
    and it is taken only when the geodesic midpoint lies on the surface.
    Resampling keeps a sample when dropping it would leave a chord through the
    mesh. Where the surface bends, the curve is sampled more finely so a
    straight segment does not cut through the interior or through empty space.
    """
    source_bary = sanitize(np.asarray(source_bary, dtype=np.float64))
    source_point = point_at(mesh, source_face, source_bary)
    bary = sanitize(np.asarray(start_bary, dtype=np.float64))
    points = [point_at(mesh, start_face, bary)]

    def finish(samples):
        return shorten_curve(mesh, samples, int(start_face))

    if np.linalg.norm(points[0] - source_point) <= 1e-8 * mesh["diagonal"]:
        return finish([points[0], source_point])

    face = int(start_face)
    carried = None
    limit = max(8, 6 * len(mesh["faces"]))
    for _ in range(limit):
        if contains_point(mesh, face, source_point):
            append_point(points, source_point)
            return finish(points)
        moved = descend_step(mesh, distance, face, bary, carried)
        if moved is None:
            break
        face, bary, carried = moved
        append_point(points, point_at(mesh, face, bary))
    if np.linalg.norm(points[-1] - source_point) <= mesh["mean_edge"]:
        append_point(points, source_point)
    return finish(points)


def curve_length(curve):
    total = 0.0
    for start, end in zip(curve, curve[1:]):
        total += float(np.linalg.norm(end - start))
    return total


def curve_hues(count):
    if count <= 0:
        return []
    return [colorsys.hsv_to_rgb(index / count, 0.75, 0.95) for index in range(count)]


def register_curves(curves):
    if ps.has_curve_network("geodesics"):
        ps.remove_curve_network("geodesics")
    if ps.has_point_cloud("geodesic samples"):
        ps.remove_point_cloud("geodesic samples")
    if not curves:
        return
    nodes = []
    edges = []
    edge_colors = []
    node_colors = []
    for curve, color in zip(curves, curve_hues(len(curves))):
        base = len(nodes)
        nodes.extend(curve)
        node_colors.extend([color] * len(curve))
        for index in range(len(curve) - 1):
            edges.append((base + index, base + index + 1))
            edge_colors.append(color)
    network = ps.register_curve_network(
        "geodesics", np.asarray(nodes), np.asarray(edges, dtype=np.int32)
    )
    network.add_color_quantity(
        "color", np.asarray(edge_colors, dtype=np.float64), defined_on="edges", enabled=True
    )
    network.set_radius(0.0015)
    # A sphere at every sample, larger than the stroke, so a vertex inside the
    # mesh is visible once the surface is hidden.
    samples = ps.register_point_cloud("geodesic samples", np.asarray(nodes, dtype=np.float64))
    samples.add_color_quantity("color", np.asarray(node_colors, dtype=np.float64), enabled=True)
    samples.set_radius(0.003)


class Session:
    def __init__(self, mesh, laplacian, surface, mesh_name):
        self.mesh = mesh
        self.laplacian = laplacian
        self.surface = surface
        self.mesh_name = mesh_name
        self.source = None
        self.distance = None
        self.curves = []
        self.press = None
        self.showing_distance = False

    def clear_source(self):
        self.source = None
        self.distance = None
        self.curves = []
        self.press = None
        if ps.has_point_cloud("source"):
            ps.remove_point_cloud("source")
        if ps.has_curve_network("geodesics"):
            ps.remove_curve_network("geodesics")
        if ps.has_point_cloud("geodesic samples"):
            ps.remove_point_cloud("geodesic samples")
        if self.showing_distance:
            self.surface.remove_quantity("geodesic distance")
            self.showing_distance = False
        print("click the mesh to place a source")

    def set_source(self, face, bary):
        self.curves = []
        if ps.has_curve_network("geodesics"):
            ps.remove_curve_network("geodesics")
        if ps.has_point_cloud("geodesic samples"):
            ps.remove_point_cloud("geodesic samples")
        self.source = (face, sanitize(np.asarray(bary, dtype=np.float64)))
        self.distance = heat_distance(self.mesh, self.laplacian, face, self.source[1])
        self.surface.add_scalar_quantity(
            "geodesic distance", self.distance, defined_on="vertices", enabled=True
        )
        self.showing_distance = True
        point = point_at(self.mesh, face, self.source[1])
        if ps.has_point_cloud("source"):
            ps.remove_point_cloud("source")
        cloud = ps.register_point_cloud("source", point.reshape(1, 3))
        cloud.set_radius(0.004)
        cloud.set_color((0.85, 0.15, 0.1))
        factor = self.mesh["time"] / self.mesh["mean_edge"] ** 2 if self.mesh["mean_edge"] else 1.0
        print(
            f"source on face {face}, distance {np.min(self.distance):.4g} .. {np.max(self.distance):.4g},"
            f" heat time {factor:.3g} edge lengths squared"
        )

    def add_curve(self, face, bary):
        source_face, source_bary = self.source
        curve = trace_geodesic(self.mesh, self.distance, source_face, source_bary, face, bary)
        if len(curve) < 2 or curve_length(curve) <= 1e-8 * self.mesh["diagonal"]:
            return
        self.curves.append(curve)
        register_curves(self.curves)
        reached = sample(self.mesh, self.distance, face, sanitize(np.asarray(bary, dtype=np.float64)))
        print(
            f"geodesic {len(self.curves)}: length {curve_length(curve):.4g},"
            f" distance at click {reached:.4g}"
        )


def pick_surface_point(session, screen_coords):
    pick = ps.pick(screen_coords=screen_coords)
    if not pick.is_hit or pick.structure_name != session.mesh_name:
        return None
    if "index" not in pick.structure_data:
        return None
    face = int(pick.structure_data["index"])
    if face < 0 or face >= len(session.mesh["faces"]):
        return None
    bary = project_onto(session.mesh, face, np.asarray(pick.position, dtype=np.float64))
    if bary is None:
        return None
    return face, bary


def callback():
    session = callback.session
    psim.TextUnformatted("Click the mesh to place the source.")
    psim.TextUnformatted("Each later click adds a geodesic back to that point.")
    psim.TextUnformatted("Drag to orbit. A click is a press without much movement.")
    if psim.Button("New source"):
        session.clear_source()
    io = psim.GetIO()
    if io.WantCaptureMouse:
        session.press = None
        return
    if io.MouseClicked[0]:
        pos = io.MousePos
        session.press = (float(pos[0]), float(pos[1]))
    if io.MouseReleased[0] and session.press is not None:
        pos = io.MousePos
        released = (float(pos[0]), float(pos[1]))
        pressed = session.press
        session.press = None
        if (released[0] - pressed[0]) ** 2 + (released[1] - pressed[1]) ** 2 <= 16.0:
            hit = pick_surface_point(session, released)
            if hit is None:
                return
            if session.source is None:
                session.set_source(*hit)
            else:
                session.add_curve(*hit)


def main():
    filename = "../../../data/Solids/torus.obj" if len(argv) < 2 else argv[1]
    manifold = hmesh.load(filename)
    hmesh.triangulate(manifold)
    mesh = build_mesh(manifold)
    laplacian = build_operators(mesh)

    ps.init()
    ps.set_ground_plane_mode("shadow_only")
    corners = [list(map(int, face)) for face in mesh["faces"]]
    surface = ps.register_surface_mesh(filename, mesh["positions"], corners)
    surface.set_selection_mode("faces_only")
    callback.session = Session(mesh, laplacian, surface, filename)
    ps.set_user_callback(callback)
    print("click the mesh to place a source")
    ps.show()


if __name__ == "__main__":
    main()
