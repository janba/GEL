//
//  Viewer.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 06/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Viewer_hpp
#define Viewer_hpp

#include <vector>
#include <GEL/GLGraphics/GLViewController.h>
#include "Manifold.h"
#include "Graph.h"

namespace PyGEL {
    using namespace GLGraphics;
    using GLManifoldViewerPtr = GLViewController*;
    using GLManifoldViewer_ptr = GLManifoldViewerPtr; // C-style alias
    
    GLManifoldViewer_ptr GLManifoldViewer_new();
    void GLManifoldViewer_event_loop(bool once);
    void GLManifoldViewer_display(GLManifoldViewer_ptr self, Manifold_ptr m, Graph_ptr g, char mode, bool smooth_shading, 
                                  const std::vector<float>& bg_color, const std::vector<double>& attrib_vec, bool reset_view, bool once);
    void GLManifoldViewer_clone_controller(GLManifoldViewer_ptr self, GLManifoldViewer_ptr other);
    void GLManifoldViewer_delete(GLManifoldViewer_ptr self);
    std::vector<double> GLManifoldViewer_get_annotation_points(GLManifoldViewer_ptr self);
    void GLManifoldViewer_set_annotation_points(GLManifoldViewer_ptr self, const std::vector<double>& data);
}

#endif /* Viewer_hpp */
