ifeq (${TARGET},debug)
CXXFLAGS 	+= -g -I/usr/X11R6/include
else
CXXFLAGS 	+= -O3 -DNDEBUG -I/usr/X11R6/include
endif 

LDFLAGS    	+= -L/usr/X11R6/lib

XLIBS      	=  -lXt -lXmu -lSM -lX11 
WIN_SYS_LIBS 	=
GLLIBS     	= -lGLU -lGL 
GLUTLIBS	= -lglut
ILLIBS		= -lILUT -lIL -lILU -lSDL -ljpeg -ltiff -lpng -lmng
NUMERICS  = -llapack  -lg2c -lblas

AR		= ar -cr
DEPFLAGS	= -MM
INSTALL		= install -m 0755

GEL_PATH    = /usr/local
