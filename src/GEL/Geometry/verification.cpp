/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "verification.h"

using namespace CGLA;
using namespace std;

namespace Geometry
{

float SqrDistance (const Vec3f& rkPoint,const Triangle& rkTri)
{
    Vec3f kDiff = rkTri.get_v0() - rkPoint;
    double fA00 = sqr_length(rkTri.get_edge(0));
    double fA01 = dot(rkTri.get_edge(0),-rkTri.get_edge(2));
    double fA11 = sqr_length(rkTri.get_edge(2));
    double fB0 = dot(kDiff,rkTri.get_edge(0));
    double fB1 = dot(kDiff,-rkTri.get_edge(2));
    double fC = sqr_length(kDiff);
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fS = fA01*fB1-fA11*fB0;
    double fT = fA01*fB0-fA00*fB1;
    double fSqrDist;

    if ( fS + fT <= fDet )
    {
        if ( fS < (double)0.0 )
        {
            if ( fT < (double)0.0 )  // region 4
            {
                if ( fB0 < (double)0.0 )
                {
                    fT = (double)0.0;
                    if ( -fB0 >= fA00 )
                    {
                        fS = (double)1.0;
                        fSqrDist = fA00+((double)2.0)*fB0+fC;
                    }
                    else
                    {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
                    }
                }
                else
                {
                    fS = (double)0.0;
                    if ( fB1 >= (double)0.0 )
                    {
                        fT = (double)0.0;
                        fSqrDist = fC;
                    }
                    else if ( -fB1 >= fA11 )
                    {
                        fT = (double)1.0;
                        fSqrDist = fA11+((double)2.0)*fB1+fC;
                    }
                    else
                    {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
                    }
                }
            }
            else  // region 3
            {
                fS = (double)0.0;
                if ( fB1 >= (double)0.0 )
                {
                    fT = (double)0.0;
                    fSqrDist = fC;
                }
                else if ( -fB1 >= fA11 )
                {
                    fT = (double)1.0;
                    fSqrDist = fA11+((double)2.0)*fB1+fC;
                }
                else
                {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                }
            }
        }
        else if ( fT < (double)0.0 )  // region 5
        {
            fT = (double)0.0;
            if ( fB0 >= (double)0.0 )
            {
                fS = (double)0.0;
                fSqrDist = fC;
            }
            else if ( -fB0 >= fA00 )
            {
                fS = (double)1.0;
                fSqrDist = fA00+((double)2.0)*fB0+fC;
            }
            else
            {
                fS = -fB0/fA00;
                fSqrDist = fB0*fS+fC;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double fInvDet = ((double)1.0)/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDist = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
        }
    }
    else
    {
        double fTmp0, fTmp1, fNumer, fDenom;

        if ( fS < (double)0.0 )  // region 2
        {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if ( fTmp1 > fTmp0 )
            {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
                {
                    fS = (double)1.0;
                    fT = (double)0.0;
                    fSqrDist = fA00+((double)2.0)*fB0+fC;
                }
                else
                {
                    fS = fNumer/fDenom;
                    fT = (double)1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
                }
            }
            else
            {
                fS = (double)0.0;
                if ( fTmp1 <= (double)0.0 )
                {
                    fT = (double)1.0;
                    fSqrDist = fA11+((double)2.0)*fB1+fC;
                }
                else if ( fB1 >= (double)0.0 )
                {
                    fT = (double)0.0;
                    fSqrDist = fC;
                }
                else
                {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                }
            }
        }
        else if ( fT < (double)0.0 )  // region 6
        {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if ( fTmp1 > fTmp0 )
            {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-((double)2.0)*fA01+fA11;
                if ( fNumer >= fDenom )
                {
                    fT = (double)1.0;
                    fS = (double)0.0;
                    fSqrDist = fA11+((double)2.0)*fB1+fC;
                }
                else
                {
                    fT = fNumer/fDenom;
                    fS = (double)1.0 - fT;
                    fSqrDist = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
                }
            }
            else
            {
                fT = (double)0.0;
                if ( fTmp1 <= (double)0.0 )
                {
                    fS = (double)1.0;
                    fSqrDist = fA00+((double)2.0)*fB0+fC;
                }
                else if ( fB0 >= (double)0.0 )
                {
                    fS = (double)0.0;
                    fSqrDist = fC;
                }
                else
                {
                    fS = -fB0/fA00;
                    fSqrDist = fB0*fS+fC;
                }
            }
        }
        else  // region 1
        {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if ( fNumer <= (double)0.0 )
            {
                fS = (double)0.0;
                fT = (double)1.0;
                fSqrDist = fA11+((double)2.0)*fB1+fC;
            }
            else
            {
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
                {
                    fS = (double)1.0;
                    fT = (double)0.0;
                    fSqrDist = fA00+((double)2.0)*fB0+fC;
                }
                else
                {
                    fS = fNumer/fDenom;
                    fT = (double)1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
                }
            }
        }
    }
    return fabs(fSqrDist);
}

}