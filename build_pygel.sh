#!/bin/sh
mkdir -p build
cd build 
cmake ..
make -j 12 
cd .. 
rm -fr dist 
python -m build -w
pip uninstall --yes PyGEL3D 
pip install dist/*whl
