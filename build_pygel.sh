#!/bin/sh
mkdir -p build
cd build 
cmake ..
make -j 12 
cd .. 
rm -fr dist 
python setup.py bdist_wheel 
pip uninstall --yes PyGEL3D 
pip install dist/*
