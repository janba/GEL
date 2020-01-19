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

using IntVector = vector<size_t>;

IntVector_ptr IntVector_new(size_t s) {
    return reinterpret_cast<IntVector_ptr>(new IntVector(s));
}

size_t IntVector_size(IntVector_ptr self) {
    return reinterpret_cast<IntVector*>(self)->size();
}


void IntVector_delete(IntVector_ptr self) {
    delete reinterpret_cast<IntVector*>(self);
}

size_t IntVector_get(IntVector_ptr self, size_t idx) {
    return (*reinterpret_cast<IntVector*>(self))[idx];
}
