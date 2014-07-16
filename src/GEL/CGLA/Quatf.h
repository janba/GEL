/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Quatf.h
 * @brief float based quaternion class
 */

#ifndef __CGLA_QUATF_H__
#define __CGLA_QUATF_H__

#include "ArithQuat.h"
#include "Vec3f.h"
#include "Vec4f.h"
#include "Mat3x3f.h"
#include "Mat4x4f.h"


namespace CGLA {

  /** \brief A float based Quaterinion class. 

	Quaternions are algebraic entities useful for rotation. */

	class Quatf : public ArithQuat<float,Vec3f,Quatf>
  {
  public:
		Quatf() : ArithQuat<float, Vec3f, Quatf>() {}

    /// Construct quaternion from vector and scalar
    Quatf(const Vec3f& imaginary, float real = 1.0f) : ArithQuat<float, Vec3f, Quatf>(imaginary, real) {}

    /// Construct quaternion from four scalars
    Quatf(float x, float y, float z, float _qw) :  ArithQuat<float, Vec3f, Quatf>(x,y,z,_qw) {}

		/// Construct quaternion from a 4D vector
		explicit Quatf(const Vec4f& v) : ArithQuat<float, Vec3f, Quatf>(v[0], v[1], v[2], v[3]) {}

	
		/// Get a 3x3 rotation matrix from a quaternion
    Mat3x3f get_Mat3x3f() const
    {
      float s = 2/norm();
      // note that the all q_*q_ are used twice (optimize)
      return Mat3x3f(Vec3f(1.0f - s*(qv[1]*qv[1] + qv[2]*qv[2]),
			         s*(qv[0]*qv[1] - qw*qv[2]),
			         s*(qv[0]*qv[2] + qw*qv[1])),
		     Vec3f(      s*(qv[0]*qv[1] + qw*qv[2]),
  			   1.0f - s*(qv[0]*qv[0] + qv[2]*qv[2]),
			         s*(qv[1]*qv[2] - qw*qv[0])),
		     Vec3f(      s*(qv[0]*qv[2] - qw*qv[1]),
			         s*(qv[1]*qv[2] + qw*qv[0]),
			   1.0f - s*(qv[0]*qv[0] + qv[1]*qv[1])));
    }

    /// Get a 4x4 rotation matrix from a quaternion
    Mat4x4f get_Mat4x4f() const
    {
      float s = 2.0f/norm();
      // note that the all q_*q_ are used twice (optimize?)
      return Mat4x4f(Vec4f(1.0f - s*(qv[1]*qv[1] + qv[2]*qv[2]),
			         s*(qv[0]*qv[1] - qw*qv[2]),
			         s*(qv[0]*qv[2] + qw*qv[1]),
		           0.0f),
		     Vec4f(      s*(qv[0]*qv[1] + qw*qv[2]),
			   1.0f - s*(qv[0]*qv[0] + qv[2]*qv[2]),
			         s*(qv[1]*qv[2] - qw*qv[0]),
			   0.0f),
		     Vec4f(      s*(qv[0]*qv[2] - qw*qv[1]),
			         s*(qv[1]*qv[2] + qw*qv[0]),
			   1.0f - s*(qv[0]*qv[0] + qv[1]*qv[1]),
			   0.0f),
		     Vec4f(0.0f, 0.0f, 0.0f, 1.0f));
    }
		
	/// Create an identity quaternion
  inline Quatf identity_Quatf()
  {
    return Quatf(Vec3f(0.0f));
  }

	};

}
#endif
