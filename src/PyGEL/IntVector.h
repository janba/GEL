//
//  IntVector.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 03/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef IntVector_hpp
#define IntVector_hpp

#include <stddef.h>
#include <stdbool.h>

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif



typedef char* IntVector_ptr;
#ifdef __cplusplus
extern "C" {
#endif
    DLLEXPORT IntVector_ptr IntVector_new(size_t s);
    DLLEXPORT size_t IntVector_get(IntVector_ptr self, size_t idx);
    DLLEXPORT size_t IntVector_size(IntVector_ptr self);
    DLLEXPORT void IntVector_delete(IntVector_ptr self);
#ifdef __cplusplus
}
#endif
#endif /* IntVector_hpp */
