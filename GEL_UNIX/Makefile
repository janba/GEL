# Very simple GNU Makefile that builds GEL on UNIX like platforms.
# It does not do dependencies but simply ploughs through all source
# files when “make” or “make GEL” is invoked on the command line.
# there is also a “clean” and a “distclean” target.
#
# Default is realease mode, but invoking make with “TARGET=debug” will
# build in debug mode. Also remember to use -j to parallelize make.
# finally note that this makefile only builds the GEL library, producing
# libGEL.a - it does not build test or apps.
#
# In the fullness of time a more complete makefile will be provided …
# or contributed! As a note, the reason for not autogenerating dependencies
# is that this requires a bit more intelligence of the makefile.
#

SHELL = /bin/sh
override SOURCEROOT		= $(shell pwd)/..

.SUFFIXES:
.SUFFIXES: .cpp .o .c

ifeq ($(strip ${CXX}),)
CXX = $(shell findcompiler.sh)
endif

OS =$(subst ${empty} ${empty},_,$(shell uname -s))
CPU =$(subst ${empty} ${empty},_,$(shell uname -m))

empty =
### The platform is determined by OS CPU and compiler 
PLATFORM = ${OS}_${CPU}_${CXX}

GEL_SUB_LIBS = CGLA HMesh Geometry GLGraphics Util

### Default target is release, unless debug is given as goal
ifndef TARGET
TARGET = release
endif

### Concatenation of platform and target yields a string used as 
### suffix of makefiles and in the name of the builddir.
PLATFORM_TARG = ${PLATFORM}_$(TARGET)

### Include platform specific makefile
include  ${SOURCEROOT}/GEL_UNIX/${PLATFORM}.mk

### BUILDDIR is the directory where the compilation for a specific 
### platform and target combination takes place.
BUILDDIR = build/$(PLATFORM_TARG)

### LIBDIR is the directory where we stick .a files pertaining to 
### the specific combination of platform and target.
LIBDIR = lib/${PLATFORM_TARG}

### BINDIR is the directory where binary files are put.
BINDIR = bin

### Append path to compiler flags
CXXFLAGS += -I${SOURCEROOT}/src/GEL/

### Generate list of source files.
sources = $(foreach lib, ${GEL_SUB_LIBS}, $(shell ls ${SOURCEROOT}/src/GEL/${lib}/*.cpp))

### Generate list if object files from list of sources.
objects = $(sources:.cpp=.o)

### Rule for building object files from C++ source files
%.o:%.cpp
	${CXX} -c ${CXXFLAGS} -o ${BUILDDIR}/$(notdir $@) $<
%.o:%.c
	${CXX} -c ${CXXFLAGS} -o ${BUILDDIR}/$(notdir $@) $<

all: GEL

clean:
	rm -fr ${BUILDDIR}

distclean:
	rm -fr build

GEL: DIR $(objects)
	${AR} ${BUILDDIR}/libGEL.a $(addprefix ${BUILDDIR}/,$(notdir $(objects)))
DIR:
	mkdir -p ${BUILDDIR}