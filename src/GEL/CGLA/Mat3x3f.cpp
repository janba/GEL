/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Mat3x3f.h"

namespace CGLA {
	Mat3x3f rotation_Mat3x3f(Axis axis, float angle)
	{
		Mat3x3f m(0.0f);

		switch(axis)
			{
			case XAXIS:
				m[0][0] = 1.0f;
				m[1][1] = cos(angle);
				m[1][2] = sin(angle);
				m[2][1] = -sin(angle);
				m[2][2] = cos(angle);
				break;
			case YAXIS:
				m[0][0] = cos(angle);
				m[0][2] = -sin(angle);
				m[2][0] = sin(angle);
				m[2][2] = cos(angle);
				m[1][1] = 1.0f;
				break;
			case ZAXIS:
				m[0][0] = cos(angle);
				m[0][1] = sin(angle);
				m[1][0] = -sin(angle);
				m[1][1] = cos(angle);
				m[2][2] = 1.0f;
				break;
			}

		return m;
	}

	Mat3x3f scaling_Mat3x3f(const Vec3f& v)
	{
		Mat3x3f m(0.0f);
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		return m;
	}


}
