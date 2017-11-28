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

using IntVector = std::vector<size_t>;

extern "C" {
    IntVector* IntVector_new(size_t s);
    size_t IntVector_get(IntVector* self, size_t size_t);
    size_t IntVector_size(IntVector* self);
    void IntVector_delete(IntVector* self);
}
#endif /* IntVector_hpp */
