/**
 * @file SimpleTrackBall.h
 * @brief Trackball based on azimuth, zenith angles.
 */
/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __GLGRAPHICS_SIMPLETRACKBALL_H__
#define __GLGRAPHICS_SIMPLETRACKBALL_H__

#include "../CGLA/CGLA.h"
#include "../CGLA/Vec3f.h"


namespace GLGraphics
{

/** \brief Simple trackball class. 

    Use it to let the mouse movement control the viewing transformation. 

		Typical usage:
		Construct the trackball as a global variable. Call the roll function
		from the GLUT mouse motion callback. Setup an idle callback which only
		calls glutPostRedisplay and call gl_view from the display function to
		set up the view transform.

		Deficiencies: This trackball has many shortcomings. For instance
		you cannot pan but only zoom and rotate. Go fix those problems!
*/
		class SimpleTrackBall
		{
				CGLA::Vec3f center;
	
				float r;     // Distance to center (origin)
				float theta; // Horizontal angle (azimuth)
				float phi;   // vertical angle (zenith)
	
				float da;    // angle increment
				float dr;    // zoom increment

				bool firsttime;   // used by the roll function
				float oldx, oldy;

				int X,Y,Z;

				void roll_x(float dx);
				void roll_y(float dy);

		public:

				/** Constructor. Call with the distance to the center. */
				SimpleTrackBall(const CGLA::Vec3f& _center, float _r): 
						center(_center),
						r(_r), theta(static_cast<float>(M_PI_2)), phi(0.0f), 
						da(static_cast<float>(M_PI/100.0)), dr(r/100.0f), firsttime(true),
						X(0), Y(1), Z(2)
						{} 

				/** Callthis to set up OpenGL viewing matrix. It will also 
						clear the view matrix. */
				void gl_view() const;

				void get_view(CGLA::Vec3f& c, CGLA::Vec3f& e, CGLA::Vec3f& u)
						{
								const CGLA::Vec3f up(0,1,0);
								float x = r * sin(theta) * cos(phi);
								float y = r * cos(theta);
								float z = r * sin(theta) * sin(phi);
								CGLA::Vec3f dir(x,y,z);
								CGLA::Vec3f eye = center + CGLA::Vec3f(dir[X], dir[Y], dir[Z]);
								c = center;
								e = eye;
								u = CGLA::Vec3f(up[X], up[Y], up[Z]);
						}

				void up_axis(char up);

				/** Move away. Typically called from the keyboard callback */
				void farther()
						{
								r += dr;
						}
		
				/** Move closer. Typically called from the keyboard callback */
				void closer()
						{
                            r = (std::max)(0.0f, r - dr);
						}

				/** Roll ball. Call with the x,y coordinates. This function is typically
						called from GLUT's mouse motion callback. */
				void roll(int x, int y);

				void set_center(const CGLA::Vec3f& _center)
						{
								center = _center;
						}
		};

}
#endif
