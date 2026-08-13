#!/bin/sh
# Build GEL/PyGEL, install the C++ package to ~/.local, then create and
# pip-install the PyGEL3D wheel.
set -eu
cd "$(dirname "$0")"

cmake -S . -B build -DCMAKE_INSTALL_PREFIX="${HOME}/.local"
cmake --build build -j 12
cmake --install build

rm -fr dist
python -m build -nwx
pip uninstall --yes PyGEL3D
pip install dist/*whl
