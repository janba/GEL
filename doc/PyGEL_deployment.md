# Making a distribution version of PyGEL #

It is not completely trivial to create a version of PyGEL for distribution. The main challenge lies in the fact that PyGEL is a Python package which depends on a C++ library which in turn depends on other libraries. 

## Compiling on MacOS ##

On MacOS, we can compile GEL and PyGEL as follows.

`cd GEL; mkdir build; cd build; cmake .. ; make -j 12`

However, note that PyGEL relies on GLFW, so that library must be installed for the above to work. Having executed these commands, we have a version of PyGEL that we can install by issuing two more commands which creates a wheel package and installs it somewhere in our Python site-package library presumably.

`python setup.py bdist_wheel`

`pip install dist/<something>.whl`

Moreover, this installation will often only work on the local machine precisely since we don’t necessarily have GLFW installed on the next machine we try.

Make sure to manually put the GLFW package inside the build folder next to `libPyGEL.dylib` and name it `libglfw.3.dylib`

Having done so, you also need to make a change to the library:

`install_name_tool -change /opt/local/lib/libglfw.3.dylib @rpath/libglfw.3.dylib libPyGEL.dylib`

This change ensures that the local version of GLFW is found so it will work also if GLFW is not installed. Going forward, it would make sense to compile GLFW statically on MacOS and have just the one dylib.

## Compiling on Windows ##

On Windows, we also use CMake. However, the resulting library will not be useful on other computers since they are very likely to have other Windows SDK’s installed. Only on another computer with the exact same Windows SDK might it work. However, if we are only interested in compiling PyGEL for our own computer that might not be a concern.

To make it useful for others, open the Visual Studio Solution file we create from the CMake file and do as follows:

1. Select GEL in solution explorer, right click and choose Options from the menu
2. From the leftmost sub-window choose 
3. configuration properties -> C/C++ -> Code Generation
4. In the right sub-window change Runtime Library to Multi-threaded (/MT)

Repeat the process for PyGEL

Finally, download the GLFW source code, and also change that to /MT as above and make sure to link against precisely the version of the library thus compiled.

Clearly, this is an arduous process, and it would be a good idea to either build the static linking into the CMake files or perhaps not use static linking but copy the relevant DLLs into the package.

## Packing it up and sending it off ##

To make the package, we simply need to put the dylib and dll files inside the build directory. If we do this on MacOS, libPyGEL.dylib should already be in place, and we just need to put the libglfw.3.dylib and the PyGEL.dll files inside the build library. Then we make the package by issuing 

`python setup.py bdist_wheel`

from the GEL root directory. To upload to PyPi, go:

`twine upload dist/PyGEL3D-<something>.whl`

and enter the login and password when prompted.

## What about Linux? ##

On Linux you need to install GLFW and follow the steps above for MacOS. It should just work. Only remember to install the GLFW-devel package and not the regular GLFW package since that is probably missing the header files.

## In the future ##

Some cross compilation in the cloud based solution would be so much nicer.