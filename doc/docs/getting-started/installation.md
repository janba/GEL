# Installation

## Install from PyPI (Recommended)

The easiest way to install PyGEL3D is using pip:

```bash
pip install PyGEL3D
```

This will install the pre-built binary package for Windows, macOS, or Linux.

### OpenGL runtime (required)

PyGEL3D's compiled library links against OpenGL (`libGL`) and the OpenGL Utility
library (`libGLU`). Pip does not install these system libraries. A typical
desktop already has them; a minimal Linux install (container, CI, server, Colab)
usually does not, and `import pygel3d` will fail until they are present.

#### Ubuntu/Debian Linux

```bash
sudo apt-get install libgl1 libglu1
```

`libgl1` and `libglu1` are virtual packages for the OpenGL and GLU runtimes.
Apt will install whatever implementation the distribution provides. On Ubuntu
that is typically the Mesa packages (`libglu1-mesa` is the same as `libglu1`).
That name refers to the package source, not to software rendering: the GPU
driver still supplies the OpenGL implementation.

To *build* GEL from source you also need the development headers:

```bash
sudo apt-get install libgl-dev libglu1-mesa-dev
```

#### macOS

OpenGL is provided by the system. No extra packages are needed.

#### Windows

OpenGL comes with the graphics-card driver. Update the driver if the viewer
fails to start.

## Google Colab

To use PyGEL3D in Google Colab, add this to your first notebook cell:

```python
!apt-get install libgl1 libglu1
!pip install PyGEL3D
```

## Building from Source

If you need to build PyGEL3D from source (for development or if pre-built binaries don't work on your system):

### Prerequisites

- CMake (version 3.25 or higher)
- A C++ compiler with C++20 support
- Python 3.11 or higher
- OpenGL and GLU development libraries (`libgl-dev` and `libglu1-mesa-dev` on Ubuntu)
- GLFW (automatically fetched by CMake)

### Clone the Repository

```bash
git clone https://github.com/janba/GEL.git
cd GEL
```

### Build with CMake

```bash
cmake -S . -B build
cmake --build build -j 8
cmake --install build
```

This installs the C++ library and headers into `~/.local`. Use
`-DCMAKE_INSTALL_PREFIX=<prefix>` to choose a different location.

### Create and Install the Python Package

```bash
python -m build -nwx
pip install dist/PyGEL3D-*.whl
```

Alternatively, use the provided build script:

```bash
sh build_install.sh
```

This script builds and installs both GEL and PyGEL.

## Verify Installation

After installation, verify that PyGEL3D is working:

```python
import pygel3d
print(pygel3d.__version__)
```

If this runs without errors, PyGEL3D is successfully installed!

## Optional Dependencies

For full functionality, you may want to install:

- **numpy**: Required for array operations (automatically installed with pip)
- **plotly**: Required for Jupyter notebook visualization
  ```bash
  pip install plotly
  ```

## Troubleshooting

### Import Errors

If you get import errors, ensure that:
1. Python can find the PyGEL3D package
2. The compiled C++ library is in the correct location
3. OpenGL libraries are installed

### OpenGL Errors

If `import pygel3d` fails with a missing `libGL` or `libGLU` on Linux, install
the runtimes (`sudo apt-get install libgl1 libglu1` on Ubuntu/Debian). If the
viewer fails to open on a desktop machine, update the graphics driver.

### Building Issues

If building from source fails:
1. Ensure all prerequisites are installed
2. Check that you have a C++20-compatible compiler
3. Try updating CMake to the latest version
4. Check the GitHub issues page for known problems
