/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../GL/glew.h"
#include "GLViewController.h"
#include "SimpleTrackBall.h"

using namespace std;
using namespace CGLA;

namespace GLGraphics
{

		void SimpleTrackBall::roll_x(float dx)
		{
				float mouse_sign = ((dx<0)?-1.0f:1.0f);
				phi += mouse_sign * da;
		}

		void SimpleTrackBall::roll_y(float dy)
		{
				float sgn = (dy<0)?1.0f:-1.0f;
				theta += sgn * da;
				theta = min(float(M_PI)-0.001f, max(0.0001f, theta));
		}

		void SimpleTrackBall::up_axis(char up)
		{
				switch(up)
				{
						case 'x':
						case 'X':
								X=1;
								Y=2;
								Z=0;
								break;
						case 'y':
						case 'Y':
								X=0;
								Y=1;
								Z=2;
								break;
						case 'z':
						case 'Z':
								X=2;
								Y=0;
								Z=1;
								break;
				}
				theta = static_cast<float>(M_PI_2);
				phi = 0;
		}

/** Call this to set up OpenGL viewing matrix. It will also 
		clear the view matrix. */
		void SimpleTrackBall::gl_view() const
		{
				const Vec3f up(0,1,0);
				float x = r * sin(theta) * cos(phi);
				float y = r * cos(theta);
				float z = r * sin(theta) * sin(phi);
				Vec3f dir(x,y,z);
				Vec3f e = center + Vec3f(dir[X], dir[Y], dir[Z]);
				glLoadIdentity();
				gluLookAt(e[0],e[1],e[2],
									center[0],center[1],center[2],
									up[X], up[Y], up[Z]);
		}

/** Roll ball. Call with the x,y coordinates. This function is typically
		called from GLUT's mouse motion callback. */
		void SimpleTrackBall::roll(int x, int y)
		{
				if(firsttime) 
						firsttime=false;
				else
				{
						float dx = x - oldx;
						float dy = y - oldy;
					
						if(dx*dx>dy*dy)
								roll_x(dx);
						else
								roll_y(dy);
					
				}
				oldx = x;
				oldy = y;
		}

}
