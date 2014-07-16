/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <algorithm>
#include "ArithVec3Float.h"
#include "Vec3f.h"
#include "Vec3d.h"

using namespace std;

namespace CGLA {

		template<class T, class V>
		void ArithVec3Float<T,V>::get_spherical(T &theta, T &phi, T &rlen ) const
		{  
				rlen = this->length();
				theta = acos((*this)[2]/rlen);    
				if ((*this)[0]>0)
						phi = atan((*this)[1]/(*this)[0]);
				else 
						if ((*this)[0]<0)
								phi = atan((*this)[1]/(*this)[0]) + M_PI;
						else 
								phi = ((*this)[1]>0) ? M_PI_2 : -1 * M_PI_2;
		}


		template<class T, class V>
		void ArithVec3Float<T,V>::set_spherical(T theta, T phi, T rlen )
		{
				(*this)[0] = rlen * sin(theta) * cos(phi);
				(*this)[1] = rlen * sin(theta) * sin(phi);
				(*this)[2] = rlen * cos(theta);
		}

		template<class T, class V>
		void orthogonal(const ArithVec3Float<T,V>& n, 
										ArithVec3Float<T,V>& b1, 
										ArithVec3Float<T,V>& b2)
		{
			onb(normalize(n), b1, b2);
		}

    template<class T, class V>
    void onb(const ArithVec3Float<T,V>& n, 
             ArithVec3Float<T,V>& b1, 
             ArithVec3Float<T,V>& b2)
    {
      if(n[2] < -0.9999999f) // Handle the singularity
      {
        b1 = V(1.0f,  0.0f, 0.0f);
        b2 = V(0.0f, -1.0f, 0.0f);
        return;
      }
      const T a = 1.0f/(1.0f + n[2]);
      const T b = n[0]*n[1]*a;
      b1 = V(b, n[1]*n[1]*a - 1.0f, n[1]);
      b2 = V(1.0f - n[0]*n[0]*a, -b, -n[0]);
    }

		template class ArithVec3Float<float, Vec3f>;
		template void orthogonal<float,Vec3f>(const ArithVec3Float<float,Vec3f>& n, 
																					ArithVec3Float<float,Vec3f>& b1, 
																					ArithVec3Float<float,Vec3f>& b2);
    template void onb<float,Vec3f>(const ArithVec3Float<float,Vec3f>& n, 
                                   ArithVec3Float<float,Vec3f>& b1, 
                                   ArithVec3Float<float,Vec3f>& b2);

		template class ArithVec3Float<double, Vec3d>;
		template void orthogonal<double,Vec3d>(const ArithVec3Float<double,Vec3d>& n,
																					 ArithVec3Float<double,Vec3d>& b1, 
																					 ArithVec3Float<double,Vec3d>& b2);
    template void onb<double,Vec3d>(const ArithVec3Float<double,Vec3d>& n, 
                                    ArithVec3Float<double,Vec3d>& b1, 
                                    ArithVec3Float<double,Vec3d>& b2);
}
