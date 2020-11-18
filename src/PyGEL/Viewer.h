//
//  Viewer.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 06/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Viewer_hpp
#define Viewer_hpp

#if defined(__APPLE__) || defined(__linux__)
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <stdbool.h>

typedef char* GLManifoldViewer_ptr;
typedef char* Manifold_ptr;
typedef char* Graph_ptr;

#ifdef __cplusplus
extern "C" {
#endif
    DLLEXPORT GLManifoldViewer_ptr GLManifoldViewer_new();
    
    DLLEXPORT void GLManifoldViewer_event_loop(bool once);
    
    DLLEXPORT void GLManifoldViewer_display(GLManifoldViewer_ptr self,
                                            Manifold_ptr m,
                                            Graph_ptr g,
                                            char mode,
                                            bool smooth_shading,
                                            float* bg_color,
                                            double* attrib_vec,
                                            bool reset_view,
                                            bool once);
    
    DLLEXPORT void GLManifoldViewer_delete(GLManifoldViewer_ptr);
    
    DLLEXPORT size_t GLManifoldViewer_get_annotation_points(GLManifoldViewer_ptr self, double** data);

    DLLEXPORT void GLManifoldViewer_set_annotation_points(GLManifoldViewer_ptr self, int n, double* data);

#ifdef __cplusplus
}
#endif

#endif /* Viewer_hpp */
