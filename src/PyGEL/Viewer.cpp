//
//  Viewer.cpp
//  PyGEL
//
//  Created by Jakob Andreas Bærentzen on 06/10/2017.
//  Copyright © 2017 Jakob Andreas Bærentzen. All rights reserved.
//

#include <map>
#include <GEL/GL/glew.h>
#include <GLFW/glfw3.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <GEL/GLGraphics/draw.h>
#include <GEL/GLGraphics/GLViewController.h>
#include "Viewer.h"

using namespace CGLA;
using namespace Geometry;
using namespace GLGraphics;
using namespace HMesh;
using namespace std;

struct DisplayParameters {
    Manifold* m_ptr = 0;
    AMGraph3D* g_ptr = 0;
    char mode;
    bool smooth_shading;
    Vec3f* bg_color;
    double* attrib_vec;
    bool reset_view;
};

class GLManifoldViewer {
    GLFWwindow* window = 0;
    std::vector<CGLA::Vec3d> annotation_points;
    bool active_annotation = false;
    bool do_pick = false;
    bool mouse_down = false;
    GLGraphics::GLViewController* glv = 0;
    GLGraphics::ManifoldRenderer* renderer = 0;
    GLuint graph_display_list = 0;
    bool escaping = false;
    float xscale, yscale;

public:
    GLManifoldViewer();
    ~GLManifoldViewer();
    
    bool was_initialized() const {return glv != 0;}
    
    DisplayParameters display_parameters;
    void display_init();
    
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


const vector<Vec3f> kelly_colors = {
    Vec3f(255, 179, 0)/255.0,
    Vec3f(128, 62, 117)/255.0,
    Vec3f(255, 104, 0)/255.0,
    Vec3f(166, 189, 215)/255.0,
    Vec3f(193, 0, 32)/255.0,
    Vec3f(206, 162, 98)/255.0,
    Vec3f(129, 112, 102)/255.0,
    Vec3f(0, 125, 52)/255.0,
    Vec3f(246, 118, 142)/255.0,
    Vec3f(0, 83, 138)/255.0,
    Vec3f(255, 122, 92)/255.0,
    Vec3f(83, 55, 122)/255.0,
    Vec3f(255, 142, 0)/255.0,
    Vec3f(179, 40, 81)/255.0,
    Vec3f(244, 200, 0)/255.0,
    Vec3f(127, 24, 13)/255.0,
    Vec3f(147, 170, 0)/255.0,
    Vec3f(89, 51, 21)/255.0,
    Vec3f(241, 58, 19)/255.0,
    Vec3f(35, 44, 22)/255.0
};

/** Map that maps from a window pointer to a GLManifoldViewer pointer. This map is primarily needed
    for event call backs where we know the window and need to match it to a viewer object. */
map<GLFWwindow*,GLManifoldViewer*> wv_map;

/** Here we create a renderer for a mesh */
ManifoldRenderer* render_factory(char mode, Manifold& m, bool smooth_shading, double *_attrib_vec) {
    ManifoldRenderer* renderer = 0;
    switch(mode) {
        case 'w':
            return new WireframeRenderer(m,smooth_shading);
        case 'g':
            renderer = new GlazedRenderer();
            renderer->compile_display_list(m, smooth_shading);
            return renderer;
        case 'i':
            renderer = new IsophoteLineRenderer();
            renderer->compile_display_list(m, smooth_shading);
            return renderer;
        case 's':
        {
            VertexID vid = *(m.vertices().begin());
            VertexAttributeVector<double> attrib_vec;
            double val0 = _attrib_vec[vid.get_index()];
            double min_val = val0;
            double max_val = val0;
            for(VertexID v : m.vertices())
            {
                double val = _attrib_vec[v.get_index()];
                attrib_vec[v] = val;
                min_val = min(min_val,val);
                max_val = max(max_val,val);
            }
            auto srenderer = new ScalarFieldRenderer();
            srenderer->compile_display_list(m, smooth_shading, attrib_vec, min_val, max_val, 2.2,1,1,1);
            renderer = srenderer;
        }
            return renderer;
        case 'l':
        {
            VertexAttributeVector<Vec3d> attrib_vec;
            for(VertexID v : m.vertices())
            {
                size_t idx = v.get_index();
                Vec3d val(_attrib_vec[3*idx],_attrib_vec[3*idx+1],_attrib_vec[3*idx+2]);
                attrib_vec[v] = val;
            }
            LineFieldRenderer* lrenderer = new LineFieldRenderer();
            lrenderer->compile_display_list(m, attrib_vec);
            renderer = lrenderer;
        }
            return renderer;
        case 'x':
            renderer = new GhostRenderer();
            renderer->compile_display_list(m, smooth_shading);
            return renderer;
        case 'n':
        default:
            renderer = new NormalRenderer();
            renderer->compile_display_list(m, smooth_shading);
            return renderer;
    }
    return renderer;
}

void draw_colored_ball(const Vec3d& p, float r, const Vec3f& col) {
    glUseProgram(0);
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    glColor3f(col[0],col[1],col[2]);
    glPushMatrix();
    glTranslated(p[0],p[1],p[2]);
    glScaled(r,r,r);
    draw_ball();
    glPopMatrix();
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(wv_map[window]->was_initialized() == false) return;

    if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE)
        wv_map[window]->set_escaping_true();
    if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
        wv_map[window]->clear_annotation();

}


void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(wv_map[window]->was_initialized() == false) return;

    int W,H;
    glfwGetWindowSize(window, &W, &H);
    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    wv_map[window]->mouse_pos = Vec2i(xscale*xpos, yscale*(H-ypos));
    wv_map[window]->roll_ball();
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(wv_map[window]->was_initialized() == false) return;

    if(action == GLFW_PRESS) {
        if(button == GLFW_MOUSE_BUTTON_LEFT )
            wv_map[window]->grab_ball(ROTATE_ACTION);
        else if(button == GLFW_MOUSE_BUTTON_RIGHT)
        {
            if(mods & GLFW_MOD_SHIFT)
                wv_map[window]->grab_ball(PAN_ACTION);
            else
                wv_map[window]->grab_ball(ZOOM_ACTION);
        }
    }
    else {
        wv_map[window]->release_ball();
        if(mods & GLFW_MOD_CONTROL)
            wv_map[window]->set_picking_true();
    }
}

void resize_callback(GLFWwindow* window, int width, int height)
{
    if(wv_map[window]->was_initialized() == false) return;
    wv_map[window]->display_init();
}

void GLManifoldViewer::display_init() {
    glfwMakeContextCurrent(window);
    glEnable(GL_DEPTH_TEST);

    const Vec3f& c = *(display_parameters.bg_color);
    glClearColor(c[0],c[1],c[2],1.0);
    Vec3d ctr;
    float rad;
    if (display_parameters.m_ptr != 0)
        bsphere(*(display_parameters.m_ptr), ctr, rad);
    else {
        double radd;
        tie(ctr,radd) = approximate_bounding_sphere(*display_parameters.g_ptr);
        rad = radd;
    }
    int W,H,WF,HF;
    glfwGetWindowSize(window, &W, &H);
    glfwGetWindowContentScale(window, &xscale, &yscale);
    glfwGetFramebufferSize(window, &WF, &HF);
    if(glv==0 || display_parameters.reset_view) {
        delete glv;
        glv = new GLViewController(W*xscale,H*yscale,Vec3f(ctr),2.0*rad);
    }
    else glv->reshape(WF, HF);
    if(renderer != 0) {
        delete renderer;
        renderer = 0;
    }
    if(display_parameters.m_ptr != 0)
        renderer = render_factory(display_parameters.mode,
                                  *(display_parameters.m_ptr),
                                  display_parameters.smooth_shading,
                                  display_parameters.attrib_vec);
    
    if (display_parameters.g_ptr != 0) {
        graph_display_list = glGenLists(1);
        glNewList(graph_display_list, GL_COMPILE);
        draw(*display_parameters.g_ptr);
        glEndList();
    }

    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetWindowSizeCallback(window, resize_callback);
    glfwPostEmptyEvent();
}

Vec3f new_color(const Vec3f& oc) {
    Vec3f c(oc[2],oc[0],oc[1]);
    c -= Vec3f(c.min_coord());
    return normalize(Vec3f(1)-c*c);
}

void GLManifoldViewer::display() {
    glfwMakeContextCurrent(window);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glv->set_gl_modelview();
    glv->reset_projection();
    
    if(display_parameters.m_ptr != 0)
        renderer->draw();
    
    if(display_parameters.g_ptr !=0)
        glCallList(graph_display_list);

    double rad = 0.01*glv->get_eye_dist();
    if(do_pick) {
        float depth;
        int W,H,WF,HF;
        glfwGetWindowSize(window, &W, &H);
        glfwGetFramebufferSize(window, &WF, &HF);
        double x = mouse_pos[0]*double(WF)/W;
        double y = mouse_pos[1]*double(HF)/H;
        if(depth_pick(x,y, depth)) {
            Vec3d picked_point = screen2world(x,y, depth);
            bool clicked_existing = false;
            for(int i = 0;i<annotation_points.size();++i)
                if(length(annotation_points[i]-picked_point)<rad) {
                    annotation_points[i] = Vec3d(std::nan(""));
                    clicked_existing = true;
                    break;
                }
            if(!clicked_existing) {
                size_t cap = size_t(-1);
                for(size_t i = 0;i<annotation_points.size();++i)
                    if(std::isnan(annotation_points[i][0])) {
                        cap = i;
                        break;
                    }
                if(cap == -1) {
                    cap = annotation_points.size();
                    annotation_points.push_back(Vec3d());
                }
                annotation_points[cap] = picked_point;
                active_annotation = true;
            }
        }
        do_pick = false;
    }
    if(active_annotation)
    {
        int i=0;
        for(auto& pp : annotation_points) {
            draw_colored_ball(pp, rad , kelly_colors[i%20]);
            i = i+1;
        }
    }
    glfwSwapBuffers(window);
}

GLManifoldViewer::GLManifoldViewer() {
    window = glfwCreateWindow(1024, 800, "PyGEL", NULL, NULL);
    wv_map[window] = this;

    if (!window) {
        glfwTerminate();
//        cout << "Terminating" << endl;
    }
    glfwMakeContextCurrent(window);
    glewInit();
}

GLManifoldViewer::~GLManifoldViewer() {
    wv_map.erase(window);
    delete renderer;
    delete glv;
    glfwDestroyWindow(window);
}

/** Creating a static instance of this class will cause GLFW to be initialized when needed i.e.
 in the function where the static instance  is declared. Only when the program terminates is GLFW also
 terminated */
class GLFWResource {
public:
    GLFWResource() {
//        cout << "Initializing GLFW" << endl;
        if(!glfwInit())
        {
            cout << "could not init GLFW ... bailing" << endl;
            exit(0);
        }
    }
     ~GLFWResource() {
//        cout << "Terminating GLFW" << endl;
        glfwTerminate();
    }
};

void GLManifoldViewer_event_loop(bool once) {
    do {
        glfwWaitEvents();
        for(auto WV : wv_map) {
            auto viewer = WV.second;
            if(viewer->was_initialized())
                viewer->display();
            if(viewer->get_escaping())
                return;
        }
    }
    while(!once);
}


// -----------------------------------------
// C API
// -----------------------------------------

GLManifoldViewer_ptr GLManifoldViewer_new() {
    static GLFWResource glfw_resource;
    return reinterpret_cast<GLManifoldViewer_ptr>(new GLManifoldViewer);
}


void GLManifoldViewer_display(GLManifoldViewer_ptr _self,
                              Manifold_ptr _m,
                              Graph_ptr _g,
                              char mode,
                              bool smooth_shading,
                              float* _bg_color,
                              double* attrib_vec,
                              bool reset_view,
                              bool once) {
    GLManifoldViewer* self = reinterpret_cast<GLManifoldViewer*>(_self);
    self->display_parameters = {
        reinterpret_cast<Manifold*>(_m),
        reinterpret_cast<AMGraph3D*>(_g),
        mode,
        smooth_shading,
        reinterpret_cast<Vec3f*>(_bg_color),
        attrib_vec,
        reset_view
    };
    self->display_init();
    GLManifoldViewer_event_loop(once);
}

void GLManifoldViewer_delete(GLManifoldViewer_ptr _self) {
    delete reinterpret_cast<GLManifoldViewer*>(_self);
}

size_t GLManifoldViewer_get_annotation_points(GLManifoldViewer_ptr _self, double** data) {
    GLManifoldViewer* self = reinterpret_cast<GLManifoldViewer*>(_self);
    auto& ap = self->get_annotation_points();
    size_t N = ap.size();
    *data = reinterpret_cast<double*>(&ap[0]);
    return N;
}

void GLManifoldViewer_set_annotation_points(GLManifoldViewer_ptr _self, int n, double* data) {
    GLManifoldViewer* self = reinterpret_cast<GLManifoldViewer*>(_self);
    vector<Vec3d> pts(n);
    for(int i = 0; i<n;++i)
        pts[i] = Vec3d(data[3*i], data[3*i+1], data[3*i+2]);
    self->set_annotation_points(pts);
}


