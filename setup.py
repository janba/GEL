import setuptools
from glob import glob
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyGEL3D",
    version="0.0.14",
    author="Andreas Bærentzen",
    author_email="janba@dtu.dk",
    description="PyGEL 3D (Python Bindings for GEL) contains tools for polygonal mesh based geometry processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://www2.compute.dtu.dk/projects/GEL/PyGEL/",
    package_dir = {'':'src/PyGEL'},
    packages = ['PyGEL3D'],
#    py_modules = ['gel','js'],
#    packages=setuptools.find_packages(),
    install_requires = ['numpy','plotly'], 
    data_files= [('share/lib',glob('build/*.dylib')+glob('build/*.so.*')+glob('build/*.dll'))],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS"
    ],)