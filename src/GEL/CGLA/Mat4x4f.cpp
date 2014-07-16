/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Mat4x4f.h"

namespace CGLA 
{

	Mat4x4f rotation_Mat4x4f(Axis axis, float angle)
	{
		Mat4x4f m(0.0f);

		switch(axis)
			{
			case XAXIS:
				m[0][0] = 1.0f;
				m[1][1] = cos(angle);
				m[1][2] = sin(angle);
				m[2][1] = -sin(angle);
				m[2][2] = cos(angle);
				m[3][3] = 1.0f;
				break;
			case YAXIS:
				m[0][0] = cos(angle);
				m[0][2] = -sin(angle);
				m[2][0] = sin(angle);
				m[2][2] = cos(angle);
				m[1][1] = 1.0f;
				m[3][3] = 1.0f;
				break;
			case ZAXIS:
				m[0][0] = cos(angle);
				m[0][1] = sin(angle);
				m[1][0] = -sin(angle);
				m[1][1] = cos(angle);
				m[2][2] = 1.0f;
				m[3][3] = 1.0f;
				break;
			}
		return m;
	}

	Mat4x4f translation_Mat4x4f(const Vec3f& v)
	{
		Mat4x4f m(0.0f);

		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
  
		m[0][3] = v[0];
		m[1][3] = v[1];
		m[2][3] = v[2];
  
		return m;
	}

	Mat4x4f scaling_Mat4x4f(const Vec3f& v)
	{
		Mat4x4f m(0.0f);

		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[3][3] = 1.0f;
   
		return m;
	}
    
    Mat4x4f perspective_Mat4x4f(float fovy, float aspect, float zNear, float zFar){
        assert(zNear > 0);
        assert(zFar > zNear);
        assert(fovy > 0);
        float top   = tan(fovy * DEGREES_TO_RADIANS/2) * zNear;
        float right = top * aspect;
        
        Mat4x4f c(0);
        c[0][0] = zNear/right;
        c[1][1] = zNear/top;
        c[2][2] = -(zFar + zNear)/(zFar - zNear);
        c[2][3] = -2.0*zFar*zNear/(zFar - zNear);
        c[3][2] = -1.0;
        c[3][3] = 0.0;
        
        return c;
    }
    
    Mat4x4f frustum_Mat4x4f(float left, float right, float bottom, float top, float nearVal, float farVal){
        assert(nearVal > 0);
        assert(farVal > nearVal);
        assert(right > left);
        assert(top > bottom);
        Mat4x4f c(0);
        c[0][0] = 2.0 * nearVal / (right - left);
        c[0][2] = (right + left) / (right - left);
        c[1][1] = 2.0 * nearVal / (top - bottom);
        c[1][2] = (top + bottom) / (top - bottom);
        c[2][2] = -(farVal + nearVal) / (farVal - nearVal);
        c[2][3] = -2.0 * farVal * nearVal / (farVal - nearVal);
        c[3][2] = -1.0;
        c[3][3] = 0.0;
        return c;
    }
    
    Mat4x4f ortho_Mat4x4f(float left, float right, float bottom, float top, float nearVal, float farVal) {
        assert(right > left);
        assert(top > bottom);
        assert(farVal > nearVal);
        Mat4x4f c(0);
        c[0][0] = 2.0/(right - left);
        c[1][1] = 2.0/(top - bottom);
        c[2][2] = 2.0/(nearVal - farVal);
        c[3][3] = 1.0;
        c[0][3] = -(right + left)/(right - left);
        c[1][3] = -(top + bottom)/(top - bottom);
        c[2][3] = -(farVal + nearVal)/(farVal - nearVal);
        return c;
    }
    
    
    Mat4x4f ortho2D_Mat4x4f(float left, float right, float bottom, float top){
        return ortho_Mat4x4f(left, right, bottom, top, -1, 1);
    }
    
    Mat4x4f lookAt_Mat4x4f(const Vec3f& eye, const Vec3f& at, const Vec3f& up){
        Vec4f n = Vec4f(normalize(eye - at),0.0);
        
        Vec4f u = Vec4f(normalize(cross(up,Vec3f(n))),0.0);
        Vec4f v = Vec4f(normalize(cross(Vec3f(n),Vec3f(u))),0.0);
        Vec4f t = Vec4f(0.0, 0.0, 0.0, 1.0);
        Mat4x4f c = Mat4x4f(u, v, n, t);
        return c * translation_Mat4x4f( -eye );
    }

}
