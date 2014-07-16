/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file UnitVector.h
 * @brief Unitvector coded with two angles - possibly not efficient.
 */

#ifndef __CGLA_UNITVECTOR_H__
#define __CGLA_UNITVECTOR_H__

#include "CGLA.h"
#include "Vec3f.h"
#include "TableTrigonometry.h"

namespace CGLA
{
  namespace TT=TableTrigonometry;

  /** \brief The UnitVector stores a unit length vector as two angles.

  A vector stored as two (fix point) angles is much smaller than
  vector stored in the usual way. On a 32 bit architecture this
  class should take up four bytes. not too bad. */
  class UnitVector
    {
      TT::Angle theta, phi;

      void encode(const Vec3f& v)
	{
	  theta = TT::t_atan(std::sqrt(CGLA::sqr(v[0])+CGLA::sqr(v[1])), v[2]);
	  phi   = TT::t_atan(v[0],v[1]);
	}
	
    public:

      /// Construct unitvector from normal vector
      explicit UnitVector(const Vec3f& v) {encode(v);}

      /// Construct default unit vector
      explicit UnitVector(): theta(0), phi(0) {}

      /// Get theta angle
      float t() const {return TT::angle2float(theta);}

      /// Get phi angle
      float f() const {return TT::angle2float(phi);}

      /// Reconstruct Vec3f from unit vector
      operator Vec3f() const
	{
	  float costf = TT::t_cos(theta);
	  return Vec3f(TT::t_cos(phi)*costf , 
		       TT::t_sin(phi)*costf, 
		       TT::t_sin(theta));
	}

      /// Test for equality.
      bool operator==(const UnitVector& u) const
	{
	  return theta == u.theta && phi == u.phi;
	}
    };


  /// Inline output operator.
  inline std::ostream& operator<<(std::ostream& os, const UnitVector& u)
    {
      os << "<v=" << u.t() << " h=" << u.f() << ">";
      return os;
    }

}
#endif
