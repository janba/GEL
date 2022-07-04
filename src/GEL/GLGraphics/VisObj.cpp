/*
 *  VisObj.cpp
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <GEL/GLGraphics/VisObj.h>

#include <GEL/GLGraphics/Console.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>
#include <GEL/HMesh/load.h>
#include <GEL/HMesh/curvature.h>

#include <GEL/CGLA/Mat3x3d.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Vec4d.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace GLGraphics;
using namespace Geometry;

namespace GLGraphics {

void VisObj::refit(const Vec3d& _bsc, double _bsr)
{
    bsphere_center = _bsc;
    bsphere_radius = _bsr;
    view_ctrl_vec[view_ctrl_id].set_centre(Vec3f(bsphere_center));
    view_ctrl_vec[view_ctrl_id].set_eye_dist(2*bsphere_radius);
}

void VisObj::refit()
{
    bsphere(mani, bsphere_center, bsphere_radius);
    view_ctrl_vec[view_ctrl_id].set_centre(Vec3f(bsphere_center));
    view_ctrl_vec[view_ctrl_id].set_eye_dist(2*bsphere_radius);
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

template<typename T>
bool VisObj::select_entity(const Vec2i& pos,
                           vector<pair<ItemID<T>, Vec3d>>& item_vec,
                           ItemID<T> invalid_id,
                           IDSet<T>& selection_set) {
    float d;
    if(depth_pick(pos[0], pos[1], d))
    {
        Vec3d c;
        float r;
        bsphere(mani, c, r);
        auto closest = invalid_id;
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

struct AttribStats {
    double mean=0, min_val=1e300, max_val=-1e300, std_dev=0;
    
    pair<double,double> get_range() const {
        return make_pair(max(min_val, mean - 3 * std_dev), min(max_val, mean + 3 * std_dev));
        
    }
};

AttribStats attribute_statistics(const Manifold& m, const VertexAttributeVector<double>& attrib) {
    AttribStats stats;
    
    if (attrib.size()>0) {
        int cnt = 0;
        for(auto v: m.vertices()) {
            if (v.get_index() < attrib.size()) {
                stats.mean += attrib[v];
                stats.min_val = min(attrib[v], stats.min_val);
                stats.max_val = max(attrib[v], stats.max_val);
                ++cnt;
            }
        }
        stats.mean /= cnt;
        
        double variance = 0;
        for(auto v: m.vertices())
            if (v.get_index() < attrib.size()) {
                variance += sqr(attrib[v]-stats.mean);
            }
        variance /= cnt;
        
        stats.std_dev = sqrt(variance);
    }
    return stats;
}

void VisObj::produce_renderer(const std::string& display_method , Console& cs, bool smooth, float gamma)
{
    delete renderer;
    
    string short_name = display_method.substr(0,3);
    
    static Console::variable<int> use_shading(1);
    static Console::variable<int> use_stripes(1);
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
        static Console::variable<string> line_direction("max");
        static Console::variable<string> method("tensors");
        static Console::variable<int> smoothing_iter(0);
        
        line_direction.reg(cs,"display.curvature_lines.direction", "");
        method.reg(cs, "display.curvature_lines.method", "");
        smoothing_iter.reg(cs, "display.curvature_lines.smoothing_iter", "");
        
        VertexAttributeVector<Mat3x3d> curvature_tensors;
        VertexAttributeVector<Vec3d> min_curv_direction;
        VertexAttributeVector<Vec3d> max_curv_direction;
        string _line_direction = line_direction;
        VertexAttributeVector<Vec3d>& lines = (_line_direction == "min") ? min_curv_direction : max_curv_direction;
        VertexAttributeVector<Vec2d> curvature;
        
        if(string(method) == "tensors")
        {
            curvature_tensors_from_edges(mani, curvature_tensors);
            
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

        smooth_vectors_on_mesh(mani, lines, smoothing_iter);

        renderer = new LineFieldRenderer();
        dynamic_cast<LineFieldRenderer*>(renderer)->compile_display_list(mani, lines);
    }
    else if(short_name == "gau"){
        static Console::variable<float> smoothing(2.0f);
        smoothing.reg(cs, "display.gaussian_curvature_renderer.smoothing", "");
        VertexAttributeVector<double> scalars(mani.allocated_vertices());
        gaussian_curvature_angle_defects(mani, scalars, smoothing);
        auto stats = attribute_statistics(mani, scalars);
        auto [min_G, max_G] = stats.get_range();
        renderer = new ScalarFieldRenderer();
        dynamic_cast<ScalarFieldRenderer*>(renderer)->compile_display_list(mani, smooth, scalars, min_G, max_G, gamma,use_stripes,color_sign,use_shading);
        
    }
    else if(short_name == "mea"){
        static Console::variable<int> smoothing(2);
        smoothing.reg(cs, "display.mean_curvature_renderer.smoothing", "");
        
        VertexAttributeVector<double> scalars(mani.allocated_vertices());
        mean_curvatures(mani, scalars, smoothing);
        
        auto stats = attribute_statistics(mani, scalars);
        auto [min_G, max_G] = stats.get_range();
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
        auto stats = attribute_statistics(mani, scalar_field);
        auto [min_G, max_G] = stats.get_range();
        renderer = new ScalarFieldRenderer();
        dynamic_cast<ScalarFieldRenderer*>(renderer)->compile_display_list(mani, smooth,scalar_field, min_G, max_G, gamma,use_stripes,color_sign,use_shading);
    }
    else if(short_name == "col")
    {
        renderer = new ColorFieldRenderer();
        dynamic_cast<ColorFieldRenderer*>(renderer)->compile_display_list(mani,
                                                                          smooth,
                                                                          color_field,
                                                                          gamma);
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
    float r = view_ctrl_vec[view_ctrl_id].get_eye_dist();
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
            circulate_face_ccw(mani, fid, std::function<void(VertexID)>([&](VertexID v){
                glVertex3dv(mani.pos(v).get());
            }) );
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

void VisObj::construct_obb_tree()
{
    build_OBBTree(mani, obb_tree);
}


void VisObj::display(const std::string& display_method , Console& cs, bool smooth, float gamma)
{
    if(create_display_list){
        create_display_list = false;
        produce_renderer(display_method, cs, smooth, gamma);
        glDeleteLists(graph_list, 1);
        graph_list = glGenLists(1);
        glNewList(graph_list, GL_COMPILE);
        if(!graph.empty())
            draw(graph);
        glEndList();
    }
    view_ctrl_vec[view_ctrl_id].set_gl_modelview();
    renderer->draw();
    if(!vertex_selection.empty() ||
       !face_selection.empty() ||
       !halfedge_selection.empty())
        draw_selection();
    glCallList(graph_list);
    if(!obb_tree.empty())
    {
        static Console::variable<int> max_level(1);
        max_level.reg(cs, "display.obb_tree.max_level", "Maximum level for OBB Tree rendering");
        draw(obb_tree,max_level);
    }
}
}
