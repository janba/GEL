# Installation

## Install from PyPI (Recommended)

The easiest way to install PyGEL3D is using pip:

```bash
pip install PyGEL3D
```

This will install the pre-built binary package for Windows, macOS, or Linux.

### Installing OpenGL Dependencies

PyGEL3D relies on OpenGL for visualization. On most systems, OpenGL is already installed, but you may need to install additional libraries.

#### Ubuntu/Debian Linux

```bash
sudo apt-get install libglu1 libgl1
```

#### macOS

OpenGL is typically pre-installed on macOS. No additional steps are usually needed.

#### Windows

OpenGL drivers are usually installed with your graphics card drivers. If you encounter issues, update your graphics drivers.

## Google Colab

To use PyGEL3D in Google Colab, add this to your first notebook cell:

```python
!apt-get install libglu1 libgl1
!pip install PyGEL3D
```

## Building from Source

If you need to build PyGEL3D from source (for development or if pre-built binaries don't work on your system):

### Prerequisites

- CMake (version 3.15 or higher)
- A C++ compiler with C++17 support
- Python 3.7 or higher
- OpenGL development libraries
- GLFW (automatically fetched by CMake)

### Clone the Repository

```bash
git clone https://github.com/janba/GEL.git
cd GEL
```

### Build with CMake

```bash
mkdir build
cd build
cmake ..
make -j 8
sudo make install
cd ..
```

### Create and Install the Python Package

```bash
python -m build -nwx
pip install dist/PyGEL3D-*.whl
```

Alternatively, use the provided build script:

```bash
sh build_pygel.sh
```

This script automates the entire build and installation process.

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

If you encounter OpenGL-related errors:
1. Update your graphics drivers
2. Check that OpenGL is properly installed
3. On Linux, ensure X11 is configured correctly

### Building Issues

If building from source fails:
1. Ensure all prerequisites are installed
2. Check that you have a C++17-compatible compiler
3. Try updating CMake to the latest version
4. Check the GitHub issues page for known problems
