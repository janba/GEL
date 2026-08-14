# Introduction to PyGEL

*J. Andreas Bærentzen, January 2018 (revised February 2021)*

This tutorial is based on the original
[Introduction to PyGEL](https://www2.compute.dtu.dk/projects/GEL/PyGEL/intro_to_PyGEL.html)
notebook. That page includes interactive 3D views exported from Jupyter.

The PyGEL Python library for geometry processing provides Python 3 bindings for
a subset of the features of the C++ [GEL library](https://github.com/janba/GEL).

PyGEL is a Python package consisting of five modules:

- `hmesh` provides `Manifold`, which is a class that represents polygonal
  meshes using the *halfedge* representation. `hmesh` also provides a slew of
  functions for manipulating polygonal meshes and the `MeshDistance` class,
  which makes it simple to compute the distance to a triangle mesh. `hmesh`
  also includes the `skeleton_to_feq` function, which makes it straightforward
  to produce a quad-only mesh from a graph.
- `graph` contains the `Graph` class, which is used for graphs: collections of
  vertices (in 3D) connected by edges. Unlike a `Manifold`, a `Graph` does not
  have to represent a surface. Associated functions include `LS_skeleton` and
  `MSLS_skeleton`, which compute a curve skeleton from a `Graph` and return the
  result as a new `Graph`. `MSLS_skeleton` is multi-scale and hence much faster.
- `gl_display` provides the `Viewer` class, which makes it simple to visualize
  meshes and graphs.
- `jupyter_display` makes it easy to use PyGEL in a Jupyter Notebook. This
  module contains a function that creates a widget for interactively
  visualizing a mesh or a graph. It is based on Plotly, and it is possible to
  export the resulting notebooks to HTML while preserving the interactive 3D
  graphics.
- `spatial` contains the `I3DTree` class, a kD-tree specialized for mapping 3D
  points to integers — typically indices. `scipy.spatial` has a more generic
  class, so this is perhaps not the most important part of PyGEL.

For installation, see the [Installation Guide](../getting-started/installation.md).
For a concise description of the API, see the
[API Reference](../api/overview.md).

## Loading and viewing a mesh

To use PyGEL you need to import the appropriate module. The code below starts by
doing that, then it loads a mesh (a `Manifold`) called `m`, creates a `Viewer`
called `viewer`, and finally displays `m` using `viewer`. The script could be
run from a Jupyter notebook, from an editor, or from an interactive Python
shell. Adjust the mesh path to a file from the GEL `data` directory.

```python
from pygel3d import hmesh, gl_display as gl

m = hmesh.load("../data/bunnygtest.obj")
viewer = gl.Viewer()
viewer.display(m)
```

Having executed the code above, you should see a window displaying a 3D model of
a bunny, and it is now possible to play with the mesh. You can rotate by
pressing the left mouse button and dragging the mouse. You can zoom using the
right mouse button, and if you hold shift the right mouse button pans instead
of zooming.

One thing that might be puzzling is that the `display` function does not
return. That is because the viewer is in the same thread of execution as the
rest of your script: it is not running in parallel. To return to the script,
press `ESC` inside the window. Doing so, you might notice that the image
freezes but the window does not go away. That is as intended. We might want to
return to the window after all — either to visualize a different mesh or simply
to look a bit more at the one we have. To make the window active again, call
`event_loop`:

```python
viewer.event_loop()
```

You will have noticed that the window with the bunny came alive again. Hit
`ESC` one more time to exit the viewer. If you really want the window to go
away, you can either wait for the entire script to terminate or explicitly
delete the object:

```python
del viewer
```

### Exploring the Viewer

The `display` function in `Viewer` has a number of parameters which were not
involved in the simple example above. In that example, `display` was called
only with the mesh `m`. Instead of `m` we could have called `display` with a
`Graph`. The `Viewer` is happy to display both types of object even at the same
time. However, it is much more full fledged as a mesh viewer than a graph
viewer.

- `mode`: a single character that determines how the mesh is visualized:
    - `'w'` — wireframe, the default, because we so often want to see the polygons.
    - `'i'` — isophote rendering. This shows curves on the surface such that
      there is the same angle between the normal and the light source. In other
      words, isophotes are curves of even intensity of shading.
    - `'g'` — glazed (this is supposed to look a bit like ceramics).
    - `'s'` — scalar field. Draw a scalar field on the surface if you provide a
      scalar value for each vertex.
    - `'l'` — line field. The same sort of visualization as scalar field, but
      you need to provide a 3D vector for each vertex.
    - `'n'` — normal shading. Fairly boring, but sometimes this is what you want.
    - `'x'` — renders the mesh transparent in a constant color. This is
      probably only useful if you want to show a mesh and a graph together, but
      then it is very useful.
- `smooth`: if `True` (the default) we use vertex normals. Otherwise, face
  normals.
- `bg_col`: background color. The default is dark grey `[0.3, 0.3, 0.3]`.
- `data`: per-vertex data for visualization. This is either a scalar or vector
  field. It is ignored unless one of those two visualization modes is selected.
  The default is `None`.
- `reset_view`: if `False`, the view is as left in the previous display call.
  If `True`, the view is reset to the default. This can be useful if you make
  changes to a mesh and then want to return to the view you had. The default
  is `False`.
- `once`: if `True` we immediately exit the event loop and return. However, the
  window stays, and if the event loop is called from this or any other viewer,
  the window will still be responsive. Default is `False`.

#### Jupyter caveat

If you are using PyGEL from within a Jupyter Notebook you may find it
intolerable to use the `gl_display.Viewer`. The Jupyter widgets for visualizing
meshes are far more useful in that setting, and some students have insisted
that the `gl_display.Viewer` is broken. It is not, but it can feel that way
when it stops the execution of the cells in your notebook.

### Visualizing in Jupyter

Sometimes we want to show a mesh in a Jupyter Notebook. For that we need a
different visualization tool, based on Plotly. Displaying a mesh with
`jupyter_display` is as easy as with `gl_display`, but the features are fewer
and a bit different.

In the following example, `jupyter_display` is used to visualize a mesh using
wireframe (default) and flat shading. `set_export_mode` must be called if you
want to export the notebook to HTML with the interactive visualization widgets
preserved. Also, if you do not call it, `display` must be the last thing you
do in a cell.

```python
from pygel3d import jupyter_display as jd
jd.set_export_mode(True)
jd.display(m, smooth=False)
```

We will often need to associate data with the mesh. `jupyter_display` supports
only scalar fields. In the example below, the vertex *x* coordinates are used
as a scalar field.

```python
pos = m.positions()
jd.display(m, data=pos[:,0], wireframe=False)
```

Unfortunately, the Plotly-based visualization is not quite as flexible as the
one based on OpenGL, so these two examples exhaust the features of
`jupyter_display`.

## Working with meshes

The point of having a mesh representation is to be able to work with the
geometry. We need to be able to visit the faces, edges, and vertices of the
mesh and make queries about both geometry and connectivity. As a first example,
assume we want to find the smallest and the biggest *x* coordinate of any
vertex in the mesh. The code below loops over all vertices and finds the
minimum and the maximum *x* coordinate.

```python
min_x_coord = 1e32
max_x_coord = -1e32

pos = m.positions()
for v in m.vertices():
    min_x_coord = min(min_x_coord, pos[v][0])
    max_x_coord = max(max_x_coord, pos[v][0])

print(min_x_coord, max_x_coord)
```

At this point, you might be asking yourself what `pos` and `v` are. `v` is a
vertex *index*: an integer we can use to refer to a vertex. For instance, if we
have all our vertex positions in an array then `v` would be used to index into
that array, referring to a specific vertex position.

In fact, `pos` is precisely an array that contains all of the vertex positions
— and it is not a copy. Any changes that you make to `pos` are directly
reflected in the mesh. To illustrate the power of that, let us randomize all
vertex positions.

```python
from random import random
from numpy import array

pos_backup = array(pos)
for v in m.vertices():
    pos[v] = [random(), random(), random()]
jd.display(m)
```

Then restore the original positions:

```python
pos[:] = pos_backup
```

`pos_backup = array(pos)` creates a new array with the contents of `pos` and
assigns it to `pos_backup`, thereby making a copy of the positions. Afterwards,
we use the slice notation `pos[:] = pos_backup` to copy back the vertex
positions. If we had just typed `pos = pos_backup` then `pos` would simply have
been another reference to the `pos_backup` array, but the actual vertex
positions would not have been changed. Coming from C++, it can be confusing
that in Python `a = b` does not copy the contents of `b` into `a` but simply
makes `a` a new name for the data in `b`.

To give just one example of iteration over the faces of a mesh:

```python
avg_area = 0
F = m.faces()
for f in F:
    avg_area += m.area(f)
avg_area /= len(F)
print("Average area : ", avg_area)
```

Above, we exploited that `m.faces()` is just an iterable container with all the
face indices. We can do the exact same thing with halfedges.

## Circulation

Very frequently, we want to visit all the vertices that are neighbors of a
given vertex, or all the halfedges that emanate from that vertex, or all of the
faces which are incident on the vertex. Such queries can be carried out using
the `circulate_vertex` function. Given a mesh `m` and a vertex `v`, we can
find the neighbors of `v` using:

```python
N = m.circulate_vertex(v, mode='v')
```

The `'v'` indicates that we want the neighbor vertices, so when the function
returns, `N` contains an iterable sequence of the neighboring vertices. If we
had given `mode='f'` or `mode='h'`, we would have obtained the incident faces
or emanating halfedges. `mode` can be omitted and defaults to `'v'`.

In the code example below, we smooth the mesh by repeatedly assigning the
average of the neighboring vertices to each vertex. This is a crude way of
performing smoothing, but it works for demonstration purposes.

```python
from numpy import zeros

for _ in range(0, 50):
    avg_pos = zeros(pos.shape)
    for v in m.vertices():
        N = m.circulate_vertex(v, 'v')
        for vn in N:
            avg_pos[v] += pos[vn]
        avg_pos[v] /= len(N)
    pos[:] = avg_pos
jd.display(m)
```

Circulation does not only work for vertices; we can also circulate around
faces. For instance,

```python
H = m.circulate_face(f, 'h')
```

produces an iterable sequence `H` containing all the halfedges that make up the
boundary of face `f`. In a similar way, we can get the vertices of `f` and the
adjacent faces.

## Other mesh operations

The `Manifold` class has many member functions that operate on individual
halfedges, faces, or vertices and which can be used to refine or simplify the
mesh. There are also several functions outside of the `Manifold` class which
can be used to manipulate the entire mesh. See the
[HMesh API](../api/hmesh.md) for a better overview.

Below is one more simple example. We load the original bunny mesh, close its
holes, triangulate it, and simplify it to 5%. We `cleanup` to remove vertices
that are no longer used, and finally show the mesh. The numbers of vertices
before and after simplification are printed.

```python
bunny = hmesh.load("../data/bunny.obj")
print("vertices before simplification :", bunny.no_allocated_vertices())
hmesh.close_holes(bunny)
hmesh.triangulate(bunny)
hmesh.quadric_simplify(bunny, 0.05)
bunny.cleanup()
print("vertices after simplification :", bunny.no_allocated_vertices())
jd.display(bunny)
```

Clearly, the only practical thing that has been achieved by the script above is
to approximately reproduce the reduced Stanford Bunny mesh that we use in the
other examples.

## Annotating meshes

There is one more feature of `gl_display.Viewer` which is arguably important.
The viewer allows users to select *annotation points* by ctrl-clicking on the
mesh. Simply spawning a viewer as shown below is sufficient to try this out.
In the example, the annotation points are then printed. If we change the
position of an annotation point, it will also move inside the viewer.

```python
pos[:] = pos_backup
viewer = gl.Viewer()
viewer.display(m)
ap = viewer.annotation_points()
for p in ap:
    print(p)
del viewer  # or you will have a stale window lying on the desktop
```

## Finding points in space

In geometry processing, we frequently have to deal with collections of
irregularly placed points in space. To facilitate queries on this type of data,
PyGEL exposes a kD-tree class. A kD-tree allows us to search for the point
closest to a given query point much faster (with asymptotic complexity
$O(\log N)$ rather than $O(N)$) than if we simply look through the point list.

In the example below, we insert all vertices from our mesh in the `I3DTree`
data structure and then use `closest_point` to locate the vertex closest to the
origin. It is less flexible than the SciPy alternative and maybe slightly
simpler to use.

```python
from pygel3d import spatial

tree = spatial.I3DTree()

for v in m.vertices():
    tree.insert(pos[v], v)
tree.build()
k, v = tree.closest_point([0, 0, 0], 1.0)
print("key = ", k, " idx = ", v)
```

## Computing distance fields

Another frequently used representation for geometry is distance fields. A
distance field is simply a function that maps a point in space to the distance
from that point to the closest point on a given surface. `MeshDistance` allows
precisely for the computation of such signed distances from arbitrary points in
space to the closest point on a `Manifold`.

```python
mdist = hmesh.MeshDistance(m)
print("Distance from origin :", mdist.signed_distance([0, 0, 0]))
```

## Skeletonization

So far we have not touched upon the `Graph` class. The reason this feature was
included was mainly to provide easy access to our skeletonization algorithm
([Bærentzen and Rotenberg](https://arxiv.org/abs/2007.03483)). Since we have
been sticking with the bunny so far, let us turn it into a skeleton as a final
trick. The skeleton of the bunny looks a bit weird in isolation, but you should
be able to guess which edges correspond to which features.

```python
from pygel3d import graph

g = graph.from_mesh(bunny)
s = graph.LS_skeleton(g)
jd.display(s)
```

`MSLS_skeleton` is the multi-scale variant and is much faster on larger graphs.

## Next steps

- [Mesh Operations](mesh-operations.md) — loading, repair, smoothing, subdivision
- [Graph Processing](graph-processing.md) — building and processing spatial graphs
- [Visualization](visualization.md) — more viewer and Jupyter display options
- [Examples](../examples.md) — complete runnable scripts
- [HMesh API](../api/hmesh.md), [Graph API](../api/graph.md)
