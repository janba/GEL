import setuptools
from glob import glob
from shutil import copyfile, copytree
from os import path,makedirs
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    """libPyGEL is a compiled extension, so the wheel must be platform-tagged."""
    def has_ext_modules(self):
        return True

with open("README.md", "r") as fh:
    long_description = fh.read()

# Create the build directory
makedirs("build/pygel3d",exist_ok=True)

# Now copy the python files to build directory
copytree("src/PyGEL/pygel3d","build/pygel3d",dirs_exist_ok=True)

# Copy the libraries to the right place.	
libs_data = []
libs = glob('build/*.dylib')+glob('build/*.so*')+glob('build/*.dll',recursive=True)
for lib_file in libs:
    _,fn = path.split(lib_file)
    dst = "build/pygel3d/"+fn
    copyfile(lib_file,dst)
    libs_data += [fn]
print("Found these libraries: ", libs_data)

setuptools.setup(
    name="PyGEL3D",
    version="0.7.0",
    author="Andreas Baerentzen",
    author_email="janba@dtu.dk",
    description="PyGEL 3D (Python Bindings for GEL) contains tools for polygonal mesh based geometry processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/janba/GEL",
    license="MIT",
    python_requires=">=3.11",
    packages = ['pygel3d'],
    package_dir = {'':'build'},
    package_data = {'pygel3d':libs_data},
    install_requires = ['numpy', 'plotly>=5.24.1,<6', 'scipy'],
    distclass=BinaryDistribution,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux"
    ],)