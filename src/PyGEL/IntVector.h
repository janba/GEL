//
//  IntVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef IntVector_hpp
#define IntVector_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <vector>

using IntVector = std::vector<std::size_t>;

extern "C" {
    DLLEXPORT IntVector* IntVector_new(std::size_t s);
    DLLEXPORT std::size_t IntVector_get(IntVector* self, std::size_t size_t);
    DLLEXPORT std::size_t IntVector_size(IntVector* self);
    DLLEXPORT void IntVector_delete(IntVector* self);
}
#endif /* IntVector_hpp */
