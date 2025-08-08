//
//  IntVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef IntVector_hpp
#define IntVector_hpp

#include <vector>
#include <cstddef>

namespace PyGEL {
    using IntVector = std::vector<size_t>;
    using IntVector_ptr = IntVector*; // C-style alias
    
    IntVector_ptr IntVector_new(size_t s);
    size_t IntVector_get(IntVector_ptr self, size_t idx);
    size_t IntVector_size(IntVector_ptr self);
    void IntVector_delete(IntVector_ptr self);
}

#endif /* IntVector_hpp */
