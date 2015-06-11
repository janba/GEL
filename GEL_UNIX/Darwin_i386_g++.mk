## Platform specific makefile

## COMPILE OPTIONS

SYSDIR = /System/Library/Frameworks/

CXXFLAGS += -I/opt/local/include
ifeq (${TARGET},debug)
CXXFLAGS 	+= -g -bind_at_load
else
CXXFLAGS 	+= -O3 -DNDEBUG -bind_at_load
endif 

LDFLAGS += -L/opt/local/lib

## LINK OPTIONS
WIN_SYS_LIBS    = -framework Cocoa -framework ApplicationServices
GLLIBS     	= -framework OpenGL
GLUTLIBS	= -framework GLUT
XLIBS      	= -L/usr/X11R6/lib -lXt -lXmu -lSM -lX11 
ILLIBS		= -framework IL
NUMERICS = -framework vecLib

AR		= libtool -o
DEPFLAGS	= -MM
INSTALL		= install -m 0755

GEL_PATH    = /usr/local
