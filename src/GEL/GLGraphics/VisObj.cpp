/*
 *  VisObj.cpp
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "VisObj.h"

#include <GLGraphics/Console.h>
#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>
#include <HMesh/load.h>
#include <HMesh/curvature.h>

#include <CGLA/Mat3x3d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace GLGraphics;

int WINX=800, WINY=800;

namespace GLGraphics {
    
    void VisObj::refit()
    {
        bsphere(mani, bsphere_center, bsphere_radius);
        
        view_ctrl.set_centre(Vec3f(bsphere_center));
        view_ctrl.set_eye_dist(2*bsphere_radius);
    }
    
    bool VisObj::reload(string _file)
    {
        if(_file != "") file = _file;
        mani.clear();
        if(!load(file, mani))
            return false;
        refit();
        return true;
    }
    
    bool VisObj::add_mesh(string file)
    {
        if(!load(file, mani))
            return false;
        return true;
    }
    
    template<typename IDType>
    bool VisObj::select_entity(const Vec2i& pos, vector<pair<IDType, Vec3d>>& item_vec,
                               IDType invalid_id,
                               set<IDType>& selection_set) {
        float d;
        if(depth_pick(pos[0], pos[1], d))
        {
            Vec3d c;
            float r;
            bsphere(mani, c, r);
            IDType closest = invalid_id;
            double min_dist = DBL_MAX;
            for(auto item : item_vec)
            {
                const Vec3d& p = item.second;
                Vec3d wp = world2screen(p);
                if(sqr_length(Vec2d(wp[0],wp[1])-Vec2d(pos))<300)
                {
                    double dist = sqr_length(screen2world(pos[0], pos[1], d)-p);
                    if(dist < min_dist)
                    {
                        min_dist = dist;
                        closest = item.first;
                    }
                }
            }
            if(closest != invalid_id)
            {
                auto info = selection_set.insert(closest);
                if(info.second == false)
                    selection_set.erase(info.first);
                post_create_display_list();
                return true;

            }
        }
        return false;
    }
    
    bool VisObj::select_vertex(const CGLA::Vec2i& pos)
    {
        vector<pair<VertexID, Vec3d>> vertex_vec;
        for(auto v: mani.vertices())
            vertex_vec.push_back(make_pair(v, mani.pos(v)));
        return select_entity(pos, vertex_vec, InvalidVertexID, vertex_selection);
    }

    bool VisObj::select_face(const CGLA::Vec2i& pos)
    {
        vector<pair<FaceID, Vec3d>> face_vec;
        for(auto f: mani.faces())
            face_vec.push_back(make_pair(f, centre(mani, f)));
        return select_entity(pos, face_vec, InvalidFaceID, face_selection);
    }

    bool VisObj::select_halfedge(const CGLA::Vec2i& pos)
    {
        vector<pair<HalfEdgeID, Vec3d>> half_edge_vec;
        for(auto h: mani.halfedges()) {
            Walker w = mani.walker(h);
            if(w.halfedge()<w.opp().halfedge())
                half_edge_vec.push_back(make_pair(h,0.5*(mani.pos(w.vertex())+mani.pos(w.opp().vertex()))));
        }
        return select_entity(pos, half_edge_vec, InvalidHalfEdgeID, halfedge_selection);
    }


    void VisObj::produce_renderer(const std::string& display_method , Console& cs, bool smooth, float gamma)
    {
        delete renderer;
        
        string short_name = display_method.substr(0,3);
        
        static Console::variable<int> use_shading(0);
        static Console::variable<int> use_stripes(0);
        static Console::variable<int> color_sign(0);
        if(short_name=="mea"||short_name=="gau"||short_name=="sca")
        {
            use_shading.reg(cs, "display.scalar_field.use_shading", "use shading for scalar field visualization");
            use_stripes.reg(cs, "display.scalar_field.use_stripes", "use stripes for scalar field visualization");
            color_sign.reg(cs, "display.scalar_field.color_sign", "color according to sign when scalar field visualizing");

        }
        
        if(short_name== "wir")
            renderer = new WireframeRenderer(mani, smooth);
        
//        else if(short_name == "har")
//            renderer = new HarmonicsRenderer(mani, harm, cs);
        
        else if(short_name == "iso") {
            renderer = new IsophoteLineRenderer();
            renderer->compile_display_list(mani,smooth);
        }
        else if(short_name == "ref") {
            renderer = new ReflectionLineRenderer();
            renderer->compile_display_list(mani,smooth);
        }
        else if(short_name == "gla") {
            renderer = new GlazedRenderer();
            dynamic_cast<GlazedRenderer*>(renderer)->compile_display_list(mani,smooth);
        }
        
        
        else if(short_name == "too") {
            renderer = new ToonRenderer();
            renderer->compile_display_list(mani,smooth);
        }
        
        
        else if(short_name == "cur"){
            static Console::variable<string> line_direction("min");
            static Console::variable<string> method("tensors");
            static Console::variable<int> smoothing_iter(1);
            
            line_direction.reg(cs,"display.curvature_lines.direction", "");
            method.reg(cs, "display.curvature_lines.method", "");
            smoothing_iter.reg(cs, "display.curvature_lines.smoothing_iter", "");
            
            VertexAttributeVector<Mat3x3d> curvature_tensors(mani.allocated_vertices());
            VertexAttributeVector<Vec3d> min_curv_direction(mani.allocated_vertices());
            VertexAttributeVector<Vec3d> max_curv_direction(mani.allocated_vertices());
            string _line_direction = line_direction;
            VertexAttributeVector<Vec3d>& lines = (_line_direction == "min") ? min_curv_direction : max_curv_direction;
            VertexAttributeVector<Vec2d> curvature(mani.allocated_vertices());
            
            if(string(method) == "tensors")
            {
                curvature_tensors_from_edges(mani, curvature_tensors);
                for(int i=0;i<smoothing_iter; ++i)
                    smooth_curvature_tensors(mani,curvature_tensors);
                
                curvature_from_tensors(mani, curvature_tensors,
                                       min_curv_direction,
                                       max_curv_direction,
                                       curvature);
            }
            else
                curvature_paraboloids(mani,
                                      min_curv_direction,
                                      max_curv_direction,
                                      curvature);
            
            renderer = new LineFieldRenderer();
            dynamic_cast<LineFieldRenderer*>(renderer)->compile_display_list(mani, lines);
        }
        else if(short_name == "gau"){
            static Console::variable<float> smoothing(2.0f);
            smoothing.reg(cs, "display.gaussian_curvature_renderer.smoothing", "");
            VertexAttributeVector<double> scalars(mani.allocated_vertices());
            gaussian_curvature_angle_defects(mani, scalars, smoothing);
            double max_G = -1e32;
            double min_G = 1e32;
            for(VertexID v: mani.vertices()) {
                max_G = max((scalars[v]), max_G);
                min_G = min((scalars[v]), min_G);
            }
            renderer = new ScalarFieldRenderer();
            dynamic_cast<ScalarFieldRenderer*>(renderer)->compile_display_list(mani, smooth, scalars, min_G, max_G, gamma,use_stripes,color_sign,use_shading);
            
        }
        else if(short_name == "mea"){
            static Console::variable<int> smoothing(2);
            smoothing.reg(cs, "display.mean_curvature_renderer.smoothing", "");
            
            VertexAttributeVector<double> scalars(mani.allocated_vertices());
            mean_curvatures(mani, scalars, smoothing);
            double max_G = -1e32;
            double min_G = 1e32;
            for(VertexID v: mani.vertices()) {
//                cout << scalars[v] << endl;
                max_G = max((scalars[v]), max_G);
                min_G = min((scalars[v]), min_G);
            }
            renderer = new ScalarFieldRenderer();
            dynamic_cast<ScalarFieldRenderer*>(renderer)->compile_display_list(mani, smooth, scalars, min_G, max_G, gamma,use_stripes,color_sign,use_shading);
        }
        else if(short_name == "amb"){
            static Console::variable<int> smoothing(2);
            smoothing.reg(cs, "display.ambient_occlusion_renderer.smoothing", "");
            
            VertexAttributeVector<double> scalars(mani.allocated_vertices());
            mean_curvatures(mani, scalars, smoothing);
            double max_G = 0;
            
            for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
                max_G = max(abs(scalars[*v]), max_G);
            
            renderer = new AmbientOcclusionRenderer();
            dynamic_cast<AmbientOcclusionRenderer*>(renderer)->compile_display_list(mani, scalars, max_G);
            
        }
        else if(short_name == "deb")
        {
            static Console::variable<float> debug_renderer_ball_radius(0.001);
            debug_renderer_ball_radius.reg(cs, "display.debug_renderer.radius","");
            renderer = new DebugRenderer;
            dynamic_cast<DebugRenderer*>(renderer)->compile_display_list(mani, smooth, debug_renderer_ball_radius);
        }
        else if(short_name == "che")
        {
            renderer = new CheckerBoardRenderer;
            renderer->compile_display_list(mani, smooth);
        }
        else if(short_name == "sca")
        {
            double max_G = scalar_field[*mani.vertices_begin()];
            double min_G = scalar_field[*mani.vertices_begin()];
            for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v) {
                max_G = max((scalar_field[*v]), max_G);
                min_G = min((scalar_field[*v]), min_G);
            }
            renderer = new ScalarFieldRenderer();
            dynamic_cast<ScalarFieldRenderer*>(renderer)->compile_display_list(mani, smooth,scalar_field, min_G, max_G, gamma,use_stripes,color_sign,use_shading);
        }
        else if(short_name == "lin")
        {
            renderer = new LineFieldRenderer();
            dynamic_cast<LineFieldRenderer*>(renderer)->compile_display_list(mani, get_line_field_attrib_vector());
        }
        else if(short_name == "ghs")
        {
            renderer = new GhostRenderer();
            dynamic_cast<GhostRenderer*>(renderer)->compile_display_list(mani,smooth);
        }
        else {
            renderer = new NormalRenderer();
            renderer->compile_display_list(mani, smooth);
        }

    
    }
    
    void VisObj::draw_selection()
    {
//        Vec3d c;
//        float r;
//        bsphere(mani, c, r);
        float r = view_ctrl.get_eye_dist();
        r *= 0.003;
        glDisable(GL_LIGHTING);
        glDepthFunc(GL_LEQUAL);
        glColor3f(1,1,0);

        for(auto vid : vertex_selection)
            if(mani.in_use(vid))
            {
                Vec3d p = mani.pos(vid);
                glPushMatrix();
                glTranslated(p[0], p[1], p[2]);
                glScalef(r, r, r);
                draw_ball();
                glPopMatrix();
            }
        glColor3f(1,.8,0);
        for(auto fid: face_selection)
        {
            if(mani.in_use(fid))
            {
                glBegin(GL_POLYGON);
                circulate_face_ccw(mani, fid, [&](VertexID v){
                    glVertex3dv(mani.pos(v).get());
                });
                glEnd();
            }
        }
        glColor3f(1,.6,0);
        for(auto hid: halfedge_selection)
        {
            glLineWidth(10);
            if(mani.in_use(hid))
            {
                Walker w = mani.walker(hid);
                glBegin(GL_LINES);
                glVertex3dv(mani.pos(w.vertex()).get());
                glVertex3dv(mani.pos(w.opp().vertex()).get());
                glEnd();
            }
            glLineWidth(1);
        }
        glDepthFunc(GL_LESS);
        glEnable(GL_LIGHTING);
    }

    
    void VisObj::display(const std::string& display_method , Console& cs, bool smooth, float gamma)
    {
        if(create_display_list){
            create_display_list = false;
            produce_renderer(display_method, cs, smooth, gamma);
        }
        view_ctrl.set_gl_modelview();
        renderer->draw();
        if(!vertex_selection.empty() ||
           !face_selection.empty() ||
           !halfedge_selection.empty())
            draw_selection();
        if(!graph.empty())
            draw(graph);
    }
}