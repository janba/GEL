from skbuild import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PyGEL3D",
    version="0.1.0",
    author="Andreas Bærentzen",
    author_email="janba@dtu.dk",
    description="PyGEL 3D (Python Bindings for GEL) contains tools for polygonal mesh based geometry processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://www2.compute.dtu.dk/projects/GEL/PyGEL/",
    packages=['PyGEL'],
    package_dir={'':'src'},
    install_requires = ['numpy','plotly'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS"
    ],
    cmake_args=['-DOpenGL_GL_PREFERENCE=GLVND', '-DBUILD_SHARED_LIBS:BOOL=ON']
)
