import setuptools
from glob import glob
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

libs = glob('build/*.dylib')+glob('build/*.so*')+glob('build/*.dll')
print(libs)
setuptools.setup(
    name="PyGEL3D",
    version="0.0.18",
    author="Andreas Bærentzen",
    author_email="janba@dtu.dk",
    description="PyGEL 3D (Python Bindings for GEL) contains tools for polygonal mesh based geometry processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://www2.compute.dtu.dk/projects/GEL/PyGEL/",
    packages = ['PyGEL3D'],
    package_dir = {'':'src/PyGEL'},
    include_package_data=True,
    package_data = {'PyGEL3D':libs},
    install_requires = ['numpy','plotly'], 
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS"
    ],)