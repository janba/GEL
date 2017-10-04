//
//  IntVector.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "IntVector.h"
#include <iostream>
using namespace std;

IntVector* IntVector_new(size_t s) {
    return new IntVector(s);
}

void IntVector_delete(IntVector* self) {
    delete self;
}

int IntVector_get(IntVector* self, int idx) {
    return (*self)[idx];
}
