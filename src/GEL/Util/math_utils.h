//
//  math_utils.h
//  GEL
//
//  Created by Jakob Andreas Bærentzen on 24/01/2025.
//  Copyright © 2025 J. Andreas Bærentzen. All rights reserved.
//

#include <cmath>

inline double smoothstep(double a, double b, double x) {
    t = min(max((x-a)/(b-a), 0.0),1.0);
    return 3.0*t*t - 2.0*t*t*t;
}
