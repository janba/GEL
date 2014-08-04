/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Mat4x4f.h
 * @brief 4x4 float matrix class
 */

#ifndef __CGLA_MAT4X4_H__
#define __CGLA_MAT4X4_H__

#include "ExceptionStandard.h"
#include "CGLA.h"
#include "Vec3f.h"
#include "Vec4f.h"
#include "ArithSqMat4x4Float.h"


namespace CGLA 
{

  /** \brief 4x4 float matrix.
      This class is useful for transformations such as perspective projections 
      or translation where 3x3 matrices do not suffice. */
  class Mat4x4f: public ArithSqMat4x4Float<Vec4f, Mat4x4f>
    {
    public:
  
      /// Construct a Mat4x4f from four Vec4f vectors
      Mat4x4f(Vec4f _a, Vec4f _b, Vec4f _c, Vec4f _d): 
	ArithSqMat4x4Float<Vec4f, Mat4x4f> (_a,_b,_c,_d) {}
  
      /// Construct the NaN matrix
      Mat4x4f() {}

      /// Construct a matrix with identical elements.
      explicit Mat4x4f(float a) : ArithSqMat4x4Float<Vec4f, Mat4x4f> (a) {}
    };

  /// Create a rotation _matrix. Rotates about one of the major axes.
  Mat4x4f rotation_Mat4x4f(CGLA::Axis axis, float angle);

  /// Create a translation matrix
  Mat4x4f translation_Mat4x4f(const Vec3f&);

  /// Create a scaling matrix.
  Mat4x4f scaling_Mat4x4f(const Vec3f&);
    
    /// Creates a perspective projection similar to gluPerspective
    /// Description from gluPerspective: perspective_Mat4x4f specifies a viewing frustum into
    /// the world coordinate system. In general, the aspect ratio in perspective_Mat4x4f
    /// should match the aspect ratio of the associated viewport. For example, aspect = 2.0
    /// means the viewer's angle of view is twice as wide in x as it is in y. If the viewport
    /// is twice as wide as it is tall, it displays the image without distortion.
    Mat4x4f perspective_Mat4x4f(float fovy, float aspect, float zNear, float zFar);
    
    /// Creates a perspective matrix similar to glFrustum
    Mat4x4f frustum_Mat4x4f(float  	left,
                            float  	right,
                            float  	bottom,
                            float  	top,
                            float  	nearVal,
                            float  	farVal);
    
    /// Creates an orthographic projection matrix (similar to glOrtho)
    Mat4x4f ortho_Mat4x4f(float left,
                          float right,
                          float bottom,
                          float top,
                          float nearVal,
                          float farVal);
    
    /// Creates a 2D orthographic projection matrix (similar to gluOrtho2D)
    Mat4x4f ortho2D_Mat4x4f(float left, float right, float bottom, float top);
    
    /// Creates a view matrix similar to gluLookAt
    Mat4x4f lookAt_Mat4x4f(const Vec3f& eye, const Vec3f& at, const Vec3f& up);

  /// Create an identity matrix.
  inline Mat4x4f identity_Mat4x4f()
    {
      return Mat4x4f(Vec4f(1.0f,0.0f,0.0f,0.0f), 
		     Vec4f(0.0f,1.0f,0.0f,0.0f), 
		     Vec4f(0.0f,0.0f,1.0f,0.0f), 
		     Vec4f(0.0f,0.0f,0.0f,1.0f));
    }

  /** Compute inverse assuming that the upper-left 3x3 sub-matrix is
      orthonormal (which is the case if the transformation is only
      a concatenation of rotations and translations).
  */
  inline Mat4x4f invert_ortho(const Mat4x4f& m)
  {
    Vec3f rx(m[0][0], m[1][0], m[2][0]);
    Vec3f ry(m[0][1], m[1][1], m[2][1]);
    Vec3f rz(m[0][2], m[1][2], m[2][2]);
    Vec3f t(m[0][3], m[1][3], m[2][3]);

    return Mat4x4f(Vec4f(rx, -dot(t, rx)),
		   Vec4f(ry, -dot(t, ry)),
		   Vec4f(rz, -dot(t, rz)),
		   Vec4f(0.0f, 0.0f, 0.0f, 1.0f));
  }
}
#endif







