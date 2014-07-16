/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file Quatd.h
 * @brief Double based quaternion class
 */

#ifndef __CGLA_QUATD_H__
#define __CGLA_QUATD_H__

#include "ArithQuat.h"
#include "Vec3d.h"
#include "Vec4d.h"
#include "Mat3x3d.h"
#include "Mat4x4d.h"


namespace CGLA {

  /** \brief A float based Quaterinion class. 

	Quaternions are algebraic entities useful for rotation. */

	class Quatd : public ArithQuat<double,Vec3d,Quatd>
  {
  public:
		Quatd() : ArithQuat<double, Vec3d, Quatd>() {}

    /// Construct quaternion from vector and scalar
    Quatd(const Vec3d& imaginary, double real = 1.0) : ArithQuat<double, Vec3d, Quatd>(imaginary, real) {}

    /// Construct quaternion from four scalars
    Quatd(double x, double y, double z, double _qw) :  ArithQuat<double, Vec3d, Quatd>(x,y,z,_qw) {}

		/// Construct quaternion from a 4D vector
		explicit Quatd(const Vec4d& v) : ArithQuat<double, Vec3d, Quatd>(v[0], v[1], v[2], v[3]) {}

	
		/// Get a 3x3 rotation matrix from a quaternion
    Mat3x3d get_Mat3x3d() const
    {
      double s = 2.0/norm();
      // note that the all q_*q_ are used twice (optimize)
      return Mat3x3d(Vec3d(1.0 - s*(qv[1]*qv[1] + qv[2]*qv[2]),
			         s*(qv[0]*qv[1] - qw*qv[2]),
			         s*(qv[0]*qv[2] + qw*qv[1])),
		     Vec3d(      s*(qv[0]*qv[1] + qw*qv[2]),
  			   1.0 - s*(qv[0]*qv[0] + qv[2]*qv[2]),
			         s*(qv[1]*qv[2] - qw*qv[0])),
		     Vec3d(      s*(qv[0]*qv[2] - qw*qv[1]),
			         s*(qv[1]*qv[2] + qw*qv[0]),
			   1.0 - s*(qv[0]*qv[0] + qv[1]*qv[1])));
    }

    /// Get a 4x4 rotation matrix from a quaternion
    Mat4x4d get_Mat4x4d() const
    {
      double s = 2/norm();
      // note that the all q_*q_ are used twice (optimize?)
      return Mat4x4d(Vec4d(1.0 - s*(qv[1]*qv[1] + qv[2]*qv[2]),
			         s*(qv[0]*qv[1] - qw*qv[2]),
			         s*(qv[0]*qv[2] + qw*qv[1]),
		           0.0),
		     Vec4d(      s*(qv[0]*qv[1] + qw*qv[2]),
			   1.0 - s*(qv[0]*qv[0] + qv[2]*qv[2]),
			         s*(qv[1]*qv[2] - qw*qv[0]),
			   0.0),
		     Vec4d(      s*(qv[0]*qv[2] - qw*qv[1]),
			         s*(qv[1]*qv[2] + qw*qv[0]),
			   1.0 - s*(qv[0]*qv[0] + qv[1]*qv[1]),
			   0.0),
		     Vec4d(0.0, 0.0, 0.0, 1.0));
    }
		
	/// Create an identity quaternion
  inline Quatd identity_Quatd()
  {
    return Quatd(Vec3d(0.0));
  }

	};

}
#endif
