/*
 *  VisObj.h
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __MESHEDIT_VISOBJ_H__
#define __MESHEDIT_VISOBJ_H__


#include <string>
#include "../GL/glew.h"
#include "../HMesh/Manifold.h"
#include "../CGLA/Vec3d.h"
#include "../GLGraphics/draw.h"
#include "../GLGraphics/Console.h"
#include "../GLGraphics/GLViewController.h"
#include "../GLGraphics/ManifoldRenderer.h"
#include "../HMesh/harmonics.h"

extern int WINX;
extern int WINY;

namespace GLGraphics {

class VisObj
{
    std::string file = "";
    GLGraphics::GLViewController view_ctrl;
    bool create_display_list = true;
    
    HMesh::Manifold mani;
    HMesh::Manifold old_mani;
    
    GLGraphics::ManifoldRenderer* renderer = nullptr;
    
    bool active_selection = false;
    HMesh::VertexAttributeVector<int> vertex_selection;
    HMesh::VertexAttributeVector<double> scalar_field;
    HMesh::VertexAttributeVector<CGLA::Vec3d> line_field;
    
    HMesh::Harmonics* harm;
    CGLA::Vec3d bsphere_center;
    float bsphere_radius;
    
    void produce_renderer(const std::string& display_method , Console& cs, bool smooth, float gamma);
    void draw_selection();
public:
    
    VisObj(): view_ctrl(WINX,WINY, CGLA::Vec3f(0), 1.0) {}
    
    HMesh::VertexAttributeVector<int>& get_vertex_selection() {
        return vertex_selection;
    }
    bool select_vertex(const CGLA::Vec2i& pos);
    void clear_selection() {
        for(auto vid : mani.vertices()) vertex_selection[vid] = 0;
        active_selection = false;
    }
    
    HMesh::VertexAttributeVector<double>& get_scalar_field_attrib_vector() {
        return scalar_field;
    }

    HMesh::VertexAttributeVector<CGLA::Vec3d>& get_line_field_attrib_vector() {
        return line_field;
    }
    
    const std::string& file_name() const {return file;}
    
    float get_bsphere_radius() const { return bsphere_radius;}
    
    HMesh::Manifold& mesh() {return mani;}
    const HMesh::Manifold& mesh_old() const {return old_mani;}
    
    void save_old() {old_mani = mani;}
    void restore_old() {mani = old_mani;}
    
    GLGraphics::GLViewController& view_control() {return view_ctrl;}
    
    void refit();
    
    bool reload(std::string _file);
    
    bool add_mesh(std::string _file);
    
    void display(const std::string& display_method , GLGraphics::Console& cs, bool smooth, float gamma);
    
    void post_create_display_list()
    {
        create_display_list = true;
    }
    
    void harmonics_analyze() {
        harm = new HMesh::Harmonics(mani);
    }

    void harmonics_reset_shape() {
        harm->reset_shape();
    }

    void harmonics_partial_reconstruct(int E0, int E1, float scale) {
        harm->partial_reconstruct(E0,E1, scale);
    }

};

}
#endif
