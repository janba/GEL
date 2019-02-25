import setuptools
from glob import glob
from shutil import copyfile
from setuptools import setup
from os import path,makedirs

with open("README.md", "r") as fh:
    long_description = fh.read()

# Create the build directory
makedirs("build/PyGEL3D",exist_ok=True)

# Now copy the python files to build directory
for py_file in glob("src/PyGEL/PyGEL3D/*.py"):
    _,fn = path.split(py_file)
    copyfile(py_file,"build/PyGEL3D/"+fn)

# Copy the libraries to the right place.	
libs_data = []
libs = glob('build/*.dylib')+glob('build/*.so*')+glob('build/*.dll')
for lib_file in libs:
    _,fn = path.split(lib_file)
    dst = "build/PyGEL3D/"+fn
    copyfile(lib_file,dst)
    libs_data += [fn]
print("Found these libraries: ", libs_data)

setuptools.setup(
    name="PyGEL3D",
    version="0.0.19",
    author="Andreas Bærentzen",
    author_email="janba@dtu.dk",
    description="PyGEL 3D (Python Bindings for GEL) contains tools for polygonal mesh based geometry processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://www2.compute.dtu.dk/projects/GEL/PyGEL/",
    packages = ['PyGEL3D'],
    package_dir = {'':'build'},
    package_data = {'PyGEL3D':libs_data},
    install_requires = ['numpy','plotly'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS"
    ],)
