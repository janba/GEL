Version 0.8.1
- Fix the Jupyter display integration so it also works in Marimo notebooks.
- Guard PyPI releases against mismatches between the pushed `v*` tag and the
  package version metadata.

Version 0.8.0
- PyPI now ships per-platform wheels instead of one `py3-none-any` archive:
  manylinux_2_28 (x86_64 and aarch64), macOS 11+ universal2, and Windows
  AMD64. Each wheel contains only the matching native library. pip selects
  the right one; one binary still covers every Python 3.11+.
- Linux wheels are built on manylinux (glibc 2.28, AlmaLinux 8), not Ubuntu,
  so they install on Debian, Fedora, RHEL 8+, Ubuntu 20.04+, Colab, and
  similar. libGLU and X11 helpers are bundled; system OpenGL (`libGL`) is
  still required.
- Wheels and the sdist are built by GitHub Actions (`wheels.yml`) and
  published on `v*` tags. No more downloading CI artifacts to assemble a
  wheel locally.
- `python -m build --wheel` / `pip install .` invoke CMake and produce a
  platform-tagged wheel. `pygel3d.experimental` is included in the package.

Version 0.7.3
- macOS library now targets macOS 11.0 (Big Sur) instead of the build
  machine's OS, so Intel and Apple Silicon Macs older than macOS 26 can
  load it.

Version 0.7.2
- Declare support for Python 3.14. The wheel remains `py3-none-any` with
  the macOS, Linux, and Windows libraries inside.

Version 0.7.1
- Restore a single `py3-none-any` wheel (as in 0.6.1) that ships
  `libPyGEL.dylib`, `libPyGEL.so`, and `PyGEL.dll`. The 0.7.0 upload was
  tagged `cp313-macosx_26_0_arm64`, so pip skipped it on almost every
  machine and installed 0.6.1 instead.

Version 0.6.1
- Only small changes to the README

Version 0.6.0
- Added type hints to Python interface
- Added Rotation System Reconstruction to GEL and PyGEL
- Vastly improved kD-Tree for m-nearest neighbor queries.
- Added complete mesh volume and area computation for HMesh
- updated C++ interface for HMesh. It was a bit random which functions were members of Manifold and which not. Now, all functions that operate on mesh elements are members. However, for backward compatibility the old functions are retained.

Version 0.7.0
- Much improved documentation: examples, an Introduction to PyGEL tutorial,
  and a CMake package so other projects can `find_package(GEL)`
- `build_install.sh` builds GEL, installs the C++ library, and produces the
  PyGEL3D wheel
- `pygel3d.__version__` is now defined; Python 3.11+ is required
- MeshDistance returns a scalar for a single query point
- `hmesh.save` raises on unsupported formats (including PLY)
- Graph-to-mesh helpers live in `hmesh` (`graph_to_cylinders`,
  `graph_to_isosurface`); the old `graph.to_mesh_*` names still work
- Random numbers used as floats go through `gel_rand01()` so UINT_MAX is
  converted exactly
- PLY int32/uint32 range checks use 32-bit limits
- Plotly dependency is `>=5.24.1,<6`