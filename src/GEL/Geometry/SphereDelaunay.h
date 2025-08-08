/*
 *  SphereDelaunay.h
 *  gtm
 *
 *  Created by J. Andreas Bærentzen on 04/12/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SPHERE_DELAUNAY_H
#define SPHERE_DELAUNAY_H

#include <vector>
#include <GEL/CGLA/Vec.h>


std::vector<CGLA::Vec3i>  SphereDelaunay(const std::vector<CGLA::Vec3d>& vertices);

#endif
