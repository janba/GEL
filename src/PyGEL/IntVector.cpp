//
//  IntVector.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include "IntVector.h"
#include <vector>
#include <iostream>
using namespace std;

namespace PyGEL {

IntVector_ptr IntVector_new(size_t s) {
    return new IntVector(s);
}

size_t IntVector_size(IntVector_ptr self) {
    return self->size();
}

void IntVector_delete(IntVector_ptr self) {
    delete self;
}

size_t IntVector_get(IntVector_ptr self, size_t idx) {
    return (*self)[idx];
}

} // namespace PyGEL
