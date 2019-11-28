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
#include "../Geometry/Graph.h"
#include "../Geometry/build_bbtree.h"
#include "../GLGraphics/draw.h"
#include "../GLGraphics/Console.h"
#include "../GLGraphics/GLViewController.h"
#include "../GLGraphics/ManifoldRenderer.h"
//#include "../HMesh/harmonics.h"

namespace GLGraphics {

class VisObj
{
    std::string file = "";
    GLGraphics::GLViewController view_ctrl;
    bool create_display_list = true;
    GLuint graph_list=0;

    
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
    
    VisObj() {}
    
    
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
    
    GLGraphics::GLViewController& view_control() {return view_ctrl;}
    
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
