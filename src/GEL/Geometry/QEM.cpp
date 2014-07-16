/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>

#include "QEM.h"
#include "../LinAlg/LapackFunc.h"

using namespace LinAlg;
using namespace CGLA;
using namespace std;

namespace Geometry
{
	Vec3d QEM::opt_pos(double sv_thresh, const CGLA::Vec3d& p0) const
	{
        CMatrix U,S,V;
        SVD(A,U,S,V);
        CMatrix Sp(3,3,0.0);
        
        double s00 = S.get(0,0);
        double limit = sv_thresh * s00;
        
        Sp.set(0,0,1.0/s00);
        Vec3d diff(0.0);
        for(int i=1;i<3;++i)
        {
            double sii = S.get(i,i);
            if(sii < limit)
            {
                Sp.set(i,i,0.0);
                Vec3d vi(V.get(0,i), V.get(1,i), V.get(2,i));
                diff += vi*dot(vi, p0);
            }
            else
                Sp.set(i,i,1.0/sii);
        }
        
        CMatrix Ap = V * Sp * U.Transposed();
        
        CVector x(3);
        x = Ap * CVector(b*0.5);
        
        return diff-Vec3d(x[0],x[1],x[2]);
	}
    
}
