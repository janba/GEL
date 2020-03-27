//
//  bounding_sphere.hpp
//  GEL
//
//  Created by Andreas Bærentzen on 26/03/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//

#ifndef bounding_sphere_hpp
#define bounding_sphere_hpp

#include <array>
#include <vector>
#include <utility>
#include "../CGLA/CGLA.h"

std::pair<CGLA::Vec2d, double> triangle_circumcircle(const std::array<CGLA::Vec2d, 3>& pts);
std::pair<CGLA::Vec3d, double> tetrahedron_circumsphere(const std::array<CGLA::Vec3d, 4>& pts);
std::pair<CGLA::Vec3d, double> bounding_sphere(const std::vector<CGLA::Vec3d>& pts);

#endif /* bounding_sphere_hpp */
