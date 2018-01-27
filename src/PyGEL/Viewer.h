//
//  Viewer.hpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 06/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#ifndef Viewer_hpp
#define Viewer_hpp

#ifdef __APPLE__
#define DLLEXPORT __attribute__ ((visibility ("default")))
#else
#define DLLEXPORT __declspec(dllexport)
#endif

#include <GEL/HMesh/Manifold.h>

class GLManifoldViewer {
    GLFWwindow* window = 0;
    std::vector<CGLA::Vec3d> annotation_points;
    bool active_annotation = false;
    bool do_pick = false;
    bool mouse_down = false;
    GLGraphics::GLViewController* glv = 0;
    GLGraphics::ManifoldRenderer* renderer = 0;
    bool escaping = false;

public:
    GLManifoldViewer();
    ~GLManifoldViewer();
    
    void display_init(HMesh::Manifold& m,
                 char mode,
                 bool smooth_shading,
                 CGLA::Vec3f* bg_color,                 
                 double* attrib_vec,
                 bool reset_view);
    
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

};

extern "C" {
    DLLEXPORT GLManifoldViewer* GLManifoldViewer_new();
    
    DLLEXPORT void GLManifoldViewer_event_loop(bool once);
    
    DLLEXPORT void GLManifoldViewer_display(GLManifoldViewer* self,
                                  HMesh::Manifold* m,
                                  char mode,
                                  bool smooth_shading,
                                  CGLA::Vec3f* bg_color,
                                  double* attrib_vec,
                                  bool reset_view,
                                  bool once);
    
    DLLEXPORT void GLManifoldViewer_delete(GLManifoldViewer*);
    
    DLLEXPORT size_t GLManifoldViewer_get_annotation_points(GLManifoldViewer* self, double** data);

}

#endif /* Viewer_hpp */
