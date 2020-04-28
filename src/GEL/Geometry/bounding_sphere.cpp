//
//  bounding_sphere.cpp
//  GEL
//
//  Created by Andreas Bærentzen on 26/03/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//

#include <algorithm>
#include <random>
#include "bounding_sphere.h"

using namespace std;
using namespace CGLA;

pair<Vec2d, double> triangle_circumcircle(const array<Vec2d, 3>& pts) {
    Vec2d p2 = pts[1]-pts[0];
    Vec2d p3 = pts[2]-pts[0];
    double a = p2[0]*p3[1]-p3[0]*p2[1];
    double bx= - ((sqr(p2[0])+sqr(p2[1]))*p3[1] - (sqr(p3[0])+sqr(p3[1]))*p2[1]);
    double by=   ((sqr(p2[0])+sqr(p2[1]))*p3[0] - (sqr(p3[0])+sqr(p3[1]))*p2[0]);
    Vec2d c(pts[0][0]-bx/(2*a),pts[0][1]-by/(2*a));
    double r =sqrt(bx*bx+by*by)/(2.0*abs(a));
    return make_pair(c, r);
}

pair<Vec3d, double> tetrahedron_circumsphere(const array<Vec3d, 4>& pts) {
    Vec4d sos;
    for(int i=0;i<4;++i)
        sos[i] = sqr_length(pts[i]);

    Mat4x4d Ma,MDx,MDy,MDz,Mc;
    for(int i=0;i<4;++i)
        {
            Ma[i][0] = pts[i][0]; Ma[i][1] = pts[i][1];  Ma[i][2] = pts[i][2];  Ma[i][3] = 1.0;
            MDx[i][0] = sos[i];   MDx[i][1] = pts[i][1]; MDx[i][2] = pts[i][2]; MDx[i][3] = 1.0;
            MDy[i][0] = sos[i];   MDy[i][1] = pts[i][0]; MDy[i][2] = pts[i][2]; MDy[i][3] = 1.0;
            MDz[i][0] = sos[i];   MDz[i][1] = pts[i][0]; MDz[i][2] = pts[i][1]; MDz[i][3] = 1.0;
            Mc[i][0] = sos[i];    Mc[i][1] = pts[i][0];  Mc[i][2] = pts[i][1];  Mc[i][3] = pts[i][2];
        }
    double a = determinant(Ma);
    double Dx =  determinant(MDx);
    double Dy = -determinant(MDy);
    double Dz =  determinant(MDz);
    double c = determinant(Mc);
    
    Vec3d x(Dx,Dy,Dz);
    x /= 2.0 * a;
    
    double r = sqrt(sqr(Dx)+sqr(Dy)+sqr(Dz)-4*a*c)/(2.0*abs(a));
    
    return make_pair(x,r);
}

pair<Vec3d, double> b_sphere(const vector<Vec3d>& pts) {
    auto N = pts.size();

    switch(N) {
        case 1:
            return make_pair(pts[0],0);
        case 2:
            return make_pair(0.5*(pts[0]+pts[1]), 0.5*length(pts[0]-pts[1]));
        case 3: {
            Vec3d a = pts[1] - pts[0];
            Vec3d b = pts[2] - pts[0];
            Vec3d n = cross(a,b);
            n.normalize();
            Vec3d X,Y;
            orthogonal(n, X, Y);
            auto [c2d,r] = triangle_circumcircle({Vec2d(0),
                Vec2d(dot(a,X), dot(a,Y)),
                Vec2d(dot(b,X), dot(b,Y))});
            Vec3d c = pts[0];
            c += X * c2d[0] + Y * c2d[1];
            return make_pair(c, r);
        }
        case 4: {
            array<Vec3d,4> _pts = {pts[0],pts[1],pts[2],pts[3]};
            return tetrahedron_circumsphere(_pts);
        }
    }
    return make_pair(Vec3d(cgla_nan()),cgla_nan());
}

pair<Vec3d, double> Welzl(vector<Vec3d> P, vector<Vec3d> R) {
    if(P.empty() || R.size() == 4)
        return b_sphere(R);
    const Vec3d p = P.back();
    P.pop_back();
    auto [c,r] = Welzl(P,R);
    if(sqr_length(p-c) <= r*r)
        return make_pair(c,r);
    R.push_back(p);
    return Welzl(P, R);
}

pair<Vec3d, double> bounding_sphere(const vector<Vec3d>& pts) {
    return Welzl(pts, {});
}

pair<Vec3d, double> approximate_bounding_sphere(const vector<Vec3d>& _pts) {
    vector<Vec3d> pts(_pts);
    if(pts.size()>1000) {
        shuffle(pts.begin(), pts.end(), default_random_engine(0));
        pts.resize(1000);
    }
    return Welzl(pts, {});
}
