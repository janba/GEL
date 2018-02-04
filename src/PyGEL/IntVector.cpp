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

size_t IntVector_size(IntVector* self) {
    return self->size();
}


void IntVector_delete(IntVector* self) {
    delete self;
}

size_t IntVector_get(IntVector* self, size_t idx) {
    return (*self)[idx];
}
