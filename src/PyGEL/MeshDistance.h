//
//  MeshDistance.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 13/12/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef MeshDistance_hpp
#define MeshDistance_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

typedef char* MeshDistance_ptr;
typedef char* Manifold_ptr;

extern "C" {
    DLLEXPORT MeshDistance_ptr MeshDistance_new(Manifold_ptr m);

    DLLEXPORT void MeshDistance_delete(MeshDistance_ptr);
    
    DLLEXPORT float MeshDistance_signed_distance(MeshDistance_ptr self,
                                                 const float* p,
                                                 float upper);

    DLLEXPORT bool MeshDistance_ray_inside_test(MeshDistance_ptr self,
                                                 const float* p,
                                                 int no_rays);

    
}

#endif /* MeshDistance_hpp */
