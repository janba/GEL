## Platform specific makefile

## COMPILE OPTIONS

SYSDIR = /System/Library/Frameworks/

CXXFLAGS += -x c++ -arch x86_64 -std=c++11 -stdlib=libc++ -Wno-extra-tokens -Wno-c99-extensions
CXXFLAGS += -I/opt/local/include
ifeq (${TARGET},debug)
CXXFLAGS 	+= -g
else
CXXFLAGS 	+= -O3 -DNDEBUG
endif 

LDFLAGS += -L/opt/local/lib

## LINK OPTIONS
WIN_SYS_LIBS    = -framework Cocoa
GLLIBS     	= -framework OpenGL
GLUTLIBS	= -framework GLUT
XLIBS      	= 
NUMERICS = -framework Accelerate

AR		= libtool -o
DEPFLAGS	= -MM
INSTALL		= install -m 0755

GEL_PATH    = /usr/local
