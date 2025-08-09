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
#include <GEL/GL/glew.h>
#include <GLFW/glfw3.h>
#include <GEL/GLGraphics/GLViewController.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include "Manifold.h"
#include "Graph.h"

struct DisplayParameters {
    HMesh::Manifold* m_ptr = 0;
    Geometry::AMGraph3D* g_ptr = 0;
    char mode = 'n';
    bool smooth_shading = true;
    CGLA::Vec3f bg_color = CGLA::Vec3f(0.3,0.3,0.3);
    std::vector<double>* attrib_vec = nullptr; // changed from reference to pointer to allow reassignment
    bool reset_view = false;
};

class GLManifoldViewer {
    GLFWwindow* window = 0;
    std::vector<CGLA::Vec3d> annotation_points;
    bool active_annotation = false;
    bool do_pick = false;
    bool mouse_down = false;
    std::shared_ptr<GLGraphics::GLViewController> glv = nullptr;
    GLGraphics::ManifoldRenderer* renderer = 0;
    GLuint graph_display_list = 0;
    bool escaping = false;
    float xscale, yscale;

public:
    GLManifoldViewer();
    ~GLManifoldViewer();
    
    bool was_initialized() const {return glv != nullptr;}
    
    DisplayParameters display_parameters;

    void display_init();

    void clone_controller(const GLManifoldViewer* other) {
        glv = other->glv;
    }
    
    void display();
    
    CGLA::Vec2i mouse_pos;

    void roll_ball() {
        if(mouse_down)
            glv->roll_ball(mouse_pos);
    }
    void grab_ball(GLGraphics::TrackBallAction tba) {
        glv->grab_ball(tba, mouse_pos);
        mouse_down = true;
    }
    void release_ball() {
        glv->release_ball();
        mouse_down = false;
    }
    
    void set_picking_true() {
        do_pick = true;
    }
    
    void set_escaping_true() {
        escaping = true;
        
    }
    
    bool get_escaping() {
        if(escaping) {
            escaping = false;
            return true;
        }
        return false;
    }
    
    void clear_annotation() {
        annotation_points.clear();
        active_annotation = false;
    }
    
    std::vector<CGLA::Vec3d>& get_annotation_points() {
        return annotation_points;
    }
    
    void set_annotation_points(const std::vector<CGLA::Vec3d>& pts) {
        active_annotation = pts.size()>0 ? true : false;
        annotation_points = pts;
    }

};

namespace PyGEL {
    using namespace GLGraphics;
    using GLManifoldViewer_ptr = GLManifoldViewer*; // C-style alias
    
    GLManifoldViewer_ptr GLManifoldViewer_new();
    void GLManifoldViewer_event_loop(bool once);
    void GLManifoldViewer_display(GLManifoldViewer_ptr _self,
                                Manifold_ptr _m,
                                Graph_ptr _g,
                                char mode,
                                bool smooth_shading,
                                const CGLA::Vec3f& bg_color, 
                                std::vector<double>& attrib_vec, 
                                bool reset_view,
                                bool once);    
    void GLManifoldViewer_clone_controller(GLManifoldViewer_ptr self, GLManifoldViewer_ptr other);
    void GLManifoldViewer_delete(GLManifoldViewer_ptr self);

    // std::vector<double> GLManifoldViewer_get_annotation_points(GLManifoldViewer_ptr self);
    // void GLManifoldViewer_set_annotation_points(GLManifoldViewer_ptr self, const std::vector<double>& data);
}

#endif /* Viewer_hpp */
