/*
 *  VisObj.h
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef MESHEDIT_VISOBJ_H
#define MESHEDIT_VISOBJ_H


#include <string>
#include <GEL/GL/glew.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/build_bbtree.h>
#include <GEL/GLGraphics/draw.h>
#include <GEL/GLGraphics/Console.h>
#include <GEL/GLGraphics/GLViewController.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
//#include <GEL/HMesh/harmonics.h>

namespace GLGraphics {

using ViewCtrlVec = std::vector<GLGraphics::GLViewController>;

inline ViewCtrlVec& get_view_ctrl_vec() {
    static ViewCtrlVec vcv;
    return vcv;
}

class VisObj
{
    ViewCtrlVec& view_ctrl_vec;

    std::string file = "";
    bool create_display_list = true;
    GLuint graph_list=0;
    int view_ctrl_id;
    
    HMesh::Manifold mani;
    HMesh::Manifold old_mani;
    
    Geometry::AMGraph3D graph;
    
    Geometry::OBBTree obb_tree;
    
    GLGraphics::ManifoldRenderer* renderer = nullptr;
    
    HMesh::VertexSet vertex_selection;
    HMesh::HalfEdgeSet halfedge_selection;
    HMesh::FaceSet face_selection;
    
    
    HMesh::VertexAttributeVector<double> scalar_field;
    HMesh::VertexAttributeVector<CGLA::Vec3d> color_field;
    HMesh::VertexAttributeVector<CGLA::Vec3d> line_field;
    
    CGLA::Vec3d bsphere_center;
    float bsphere_radius;
    
    void produce_renderer(const std::string& display_method , Console& cs, bool smooth, float gamma);
    void draw_selection();
    
    template<typename T>
    bool select_entity(const CGLA::Vec2i& pos,
                       std::vector<std::pair<HMesh::ItemID<T>, CGLA::Vec3d>>& item_vec,
                       HMesh::ItemID<T> invalid_id,
                       HMesh::IDSet<T>& selection_set);

    
public:
    
    VisObj(): view_ctrl_vec(get_view_ctrl_vec()) {
        view_ctrl_id = view_ctrl_vec.size();
        view_ctrl_vec.push_back(GLViewController());
    }
    
    
    HMesh::VertexSet& get_vertex_selection() {
        return vertex_selection;
    }
    bool select_vertex(const CGLA::Vec2i& pos);
    void clear_vertex_selection() {
        vertex_selection.clear();
    }

    HMesh::FaceSet& get_face_selection() {
        return face_selection;
    }
    bool select_face(const CGLA::Vec2i& pos);
    void clear_face_selection() {
        face_selection.clear();
    }
    

    HMesh::HalfEdgeSet& get_halfedge_selection() {
        return halfedge_selection;
    }
    bool select_halfedge(const CGLA::Vec2i& pos);
    void clear_halfedge_selection() {
       halfedge_selection.clear();
    }

    
    
    HMesh::VertexAttributeVector<double>& get_scalar_field_attrib_vector() {
        return scalar_field;
    }
    
    HMesh::VertexAttributeVector<CGLA::Vec3d>& get_color_field_attrib_vector() {
        return color_field;
    }

    HMesh::VertexAttributeVector<CGLA::Vec3d>& get_line_field_attrib_vector() {
        return line_field;
    }
    
    std::string& get_file_name() {return file;}
    
    float get_bsphere_radius() const { return bsphere_radius;}
    
    HMesh::Manifold& mesh() {return mani;}
    const HMesh::Manifold& mesh_old() const {return old_mani;}
    
    Geometry::AMGraph3D& get_graph() {return graph;}
    
    void construct_obb_tree();
    
    void save_old() {old_mani = mani;}
    void restore_old() {mani = old_mani;}
    
    GLGraphics::GLViewController& view_control() {return view_ctrl_vec[view_ctrl_id];}
    
    void sync_view_control(const VisObj& vo) {
        view_ctrl_id = vo.view_ctrl_id;
    }
    
    void refit(const CGLA::Vec3d& _bsc, double _bsr);
    void refit();

    bool reload(std::string _file);
    
    bool add_mesh(std::string _file);
    
    void display(const std::string& display_method , GLGraphics::Console& cs, bool smooth, float gamma);
    
    void post_create_display_list()
    {
        create_display_list = true;
    }
 
};

}
#endif
