//
//  MeshEditor.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/10/13.
//
//
#include <thread>
#if !defined (WIN32)
#include <stdarg.h>
#include <unistd.h>
#endif

#include <regex>

#include <GEL/GL/glew.h>
#include <functional>
#include <GEL/GLGraphics/MeshEditor.h>
#include <string>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

#include <GEL/GLGraphics/Console.h>
#include <GEL/GLGraphics/glsl_shader.h>
#include <GEL/GLGraphics/ShadowBuffer.h>

#include <GEL/CGLA/CGLA.h>

#include <GEL/HMesh/HMesh.h>

#include <GEL/Util/Timer.h>

#include <GEL/GLGraphics/VisObj.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Util;

namespace GLGraphics {
    
    namespace {
        
        bool wantshelp(const std::vector<std::string> & args)
        {
            if(args.size() == 0)
                return false;
            
            string str = args[0];
            
            if(str=="help" || str=="HELP" || str=="Help" || str=="?")
                return true;
            
            return false;
        }

        static map<unsigned short,string> hotkey_to_string = {
            {'w', "wireframe render mode"},
            {'n', "normal render mode"},
            {'i', "isophotes render mode"},
            {'r', "reflection render mode"},
            {'t', "toon render mode"},
            {'g', "glazed render mode"},
            {'x', "ghst render mode"},
            {'a', "ambient occlusion render mode"},
            {'c', "color field render mode"},
            {'s', "scalar field render mode"},
            {'d', "debug render mode"},
            {'C', "curvature lines render mode"},
            {'M', "mean curvature render mode"},
            {'G', "gaussian curvature render mode"},
//            {27,  "toggle console"},
            {'1', "switch to mesh 1"},
            {'2', "switch to mesh 2"},
            {'3', "switch to mesh 3"},
            {'4', "switch to mesh 4"},
            {'5', "switch to mesh 5"},
            {'6', "switch to mesh 6"},
            {'7', "switch to mesh 7"},
            {'8', "switch to mesh 8"},
            {'9', "switch to mesh 9"},
            {' ', "clear selections"},
            {'f', "toggle smooth/flat shading"},
            {'R', "repeat last command"}
        };
        
        void console_list_controls(MeshEditor* me, const std::vector<std::string> & args)
        {
            me->printf("");
            me->printf("== MeshEdit Controls ==");
            me->printf("");
            me->printf("-- Mouse movement --");
            me->printf("Left button down                           : rotate");
            me->printf("Left button down + SHIFT                   : edit");
            me->printf("Middle button down (or right button + ALT) : zoom");
            me->printf("Right button down                          : pan");
            me->printf("-- Mouse button press --");
            me->printf("Left mouse click + CTRL                 : select");
            me->printf("-- Hotkeys (note some are upper case) --");
            me->printf("  ESC : toggles console");
            me->printf("  <-  : switch to mesh with lower number ");
            me->printf("  ->  : switch to mesh with higher number ");
            for(auto& h: hotkey_to_string) {
                me->printf("  '%c' : %s", h.first, h.second.c_str());
            }
        }
        
        /// Function that aligns two meshes.
        void console_align(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: align_with <src>");
                me->printf("This function aligns current mesh with src");
                me->printf("In practice the GLViewController of src is copied to current mesh.");
                me->printf("The argument is mandatory and must correspond to a valide mesh entry.");
                me->printf("Note that results might be unexpected if the meshes are not on the same scale");
            }
            
            if(args.size()>0){
                int src = console_arg(args, 0, 1);
                if(src <0 || src>= me->get_no_meshes())
                {
                    me->printf("src mesh out of range");
                    return;
                }
                me->align(src,me->get_active_no());
            }
            else
                me->printf("You must enter the mesh number that you want to align with");
        }
        
        void console_merge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: merge_with <src>");
                me->printf("merges src into current mesh");
                return;
            }
            int src = console_arg(args, 0, 2);
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            Manifold& m_src = me->get_mesh(src);
            m.merge(m_src);
            return;
        }

        void console_clear(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: clear");
                me->printf("clears current mesh");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            m.clear();
            return;
        }

        void console_quit(MeshEditor* me, const std::vector<std::string> & args)
        {
            exit(0);
        }

        void console_ridge_lines(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: ridge_lines");
                return;
            }
            
            me->save_active_mesh();
            
            Manifold& mani = me->active_mesh();
            
            VertexAttributeVector<Mat3x3d> curvature_tensors;
            VertexAttributeVector<Vec3d> min_curv_direction;
            VertexAttributeVector<Vec3d> max_curv_direction;
            VertexAttributeVector<Vec2d> curvature;
            
            curvature_paraboloids(mani,
                                  min_curv_direction,
                                  max_curv_direction,
                                  curvature);
            
            for(auto vid : mani.vertices())
            {
                Vec3d max_curv_dir = normalize(max_curv_direction[vid]);
                Vec3d min_curv_dir = normalize(min_curv_direction[vid]);
                double vid_min_pc = curvature[vid][0];
                double vid_max_pc = curvature[vid][1];
                bool ridge = true;
                bool ravine = true;
                Walker w = mani.walker(vid);
                Vec3d r(0);
                for(; !w.full_circle();w = w.circulate_vertex_ccw())
                {
                    Vec3d e = (mani.pos(w.vertex()) - mani.pos(vid));
                    
                    if(abs(dot(min_curv_dir,e)) > abs(dot(max_curv_dir,e)))
                    {
                        if(curvature[w.vertex()][0]<vid_min_pc+20)
                            ravine = false;
                        
                    }
                    else
                    {
                        if(curvature[w.vertex()][1]>vid_max_pc-20)
                            ridge = false;
                    }
                }
                DebugRenderer::vertex_colors[vid] = Vec3f(ridge,ravine,0.0);
            }
            for(auto fid : mani.faces())
                DebugRenderer::face_colors[fid] = Vec3f(.3,.3,.6);
            for(auto hid : mani.halfedges()) {
                
                Walker w = mani.walker(hid);
                Vec3f c0 = DebugRenderer::vertex_colors[w.opp().vertex()];
                Vec3f c1 = DebugRenderer::vertex_colors[w.vertex()];
                
                DebugRenderer::edge_colors[hid] = (c0==c1) ? c0 : Vec3f(0.1,0.1,0.3);
                
            }
        }

        void transform_mesh(Manifold& mani, const Mat4x4d& m)
        {
            for(VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
                mani.pos(*vid) = m.mul_3D_point(mani.pos(*vid));
        }
        
        void console_scale(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.scale sx sy sz");
                me->printf("Note: If only sx is provided, uniform scaling is applied");
                return;
            }
            
            Vec3d s;
            if(args.size() > 0){
                istringstream a0(args[0]);
                double scale;
                a0 >> scale;
                s = Vec3d(scale);
            }
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> s[1];
            }
            if(args.size() > 2){
                istringstream a0(args[2]);
                a0 >> s[2];
            }
            
            me->save_active_mesh();
            transform_mesh(me->active_mesh(),scaling_Mat4x4d(s));
        }

        void console_rotate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.rotate axis_x axis_y axis_z angle");
                return;
            }
            
            Vec3d a(1,0,0);
            double angle = 0;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> a[0];
            }
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> a[1];
            }
            if(args.size() > 2){
                istringstream a0(args[2]);
                a0 >> a[2];
            }
            if(args.size() > 3){
                istringstream a0(args[3]);
                a0 >> angle;
            }
            
            me->save_active_mesh();
            Quatd q;
            q.make_rot(angle, a);
            transform_mesh(me->active_mesh(),q.get_Mat4x4d());
        }

        
        void console_translate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.translate tx ty tz");
                me->printf("Note: recenters if no arguments are provided.");
                return;
            }
            
            Vec3d t;
            
            if(args.size()==0)
            {
                float rad;
                bsphere(me->active_mesh(), t, rad);
                t = -t;
            }
            else {
                if(args.size() > 0){
                    istringstream a0(args[0]);
                    a0 >> t[0];
                }
                if(args.size() > 1){
                    istringstream a0(args[1]);
                    a0 >> t[1];
                }
                if(args.size() > 2){
                    istringstream a0(args[2]);
                    a0 >> t[2];
                }
            }
            
            me->save_active_mesh();
            transform_mesh(me->active_mesh(),translation_Mat4x4d(t));
        }

        
        void console_merge_1_ring(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.merge_1_ring");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_vertex_selection();
            for(auto v: sel)
                if(m.in_use(v))
                    m.merge_one_ring(v);
        }
        
        
        void console_flip_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.flip_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h) && precond_flip_edge(m, h))
                    m.flip_edge(h);
        }

        void console_collapse_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.collapse_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h) && precond_collapse_edge(m, h))
                    m.collapse_edge(h,true);
        }

        void console_dissolve_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.dissolve_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h))
                    m.merge_faces(m.walker(h).face(), h);
        }
        
        void console_split_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.split_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h))
                    m.split_edge(h);
        }
        
        void console_stellate_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.stellate_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            for(auto f: sel)
                if(m.in_use(f))
                    m.split_face_by_vertex(f);
        }
        
        
        void console_triangulate_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.triangulate_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            for(auto f: sel)
                if(m.in_use(f))
                    triangulate(m, f);
        }
        
        void console_bridge_faces(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.bridge_faces");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            if(sel.size() != 2) {
                me->printf("You must select exactly two faces");
                return;
            }

            int i=0;
            FaceID f0, f1;
            for(auto f: sel)
            {
                int n = no_edges(m, f);
                if(i==0) {
                    f0 = f;
                    i+= n;
                }
                else {
                    i-= n;
                    f1 = f;
                }
            }
            if(i!= 0){
                me->printf("The selected faces must have same number of edges");
                return;
            }
            vector<VertexID> loop0;
			circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
                loop0.push_back(v);
            }) );
            
            vector<VertexID> loop1;
			circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
                loop1.push_back(v);
            }) );
            
            vector<pair<VertexID, VertexID> > connections;
            
            size_t L0= loop0.size();
            size_t L1= loop1.size();
            
            assert(L0==L1);
            
            size_t L = L0;
            
            float min_len = FLT_MAX;
            int j_off_min_len = -1;
            for(int j_off = 0; j_off < L; ++j_off)
            {
                float len = 0;
                for(int i=0;i<L;++i)
                    len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
                if(len < min_len)
                {
                    j_off_min_len = j_off;
                    min_len = len;
                }
            }
            for(int i=0;i<L;++i)
                connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
            // Merge the two one rings producing two faces.
            
            // Bridge the just created faces.
            vector<HalfEdgeID> newhalfedges = m.bridge_faces(f0, f1, connections);
        }


        void console_delete(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.merge_1_ring");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto vsel = me->get_vertex_selection();
            auto hsel = me->get_halfedge_selection();
            auto fsel = me->get_face_selection();
            for(auto v: vsel)
                if(m.in_use(v))
                    m.remove_vertex(v);
            for(auto h: hsel)
                if(m.in_use(h))
                    m.remove_edge(h);

            for(auto f: fsel)
                if(m.in_use(f))
                    m.remove_face(f);
        }
        
        void console_split_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.split_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto vsel = me->get_vertex_selection();
            auto fsel = me->get_face_selection();
            for(auto f: fsel)
                if(m.in_use(f))
                {
                    VertexID v0=InvalidVertexID,v1=InvalidVertexID;
					circulate_face_ccw(m, f, std::function<void(VertexID)>( [&](VertexID v){
                        if(vsel.count(v))
                        {
                            if(v0==InvalidVertexID)
                                v0 = v;
                            else
                                v1 = v;
                        
                        }
                    
                    }) );
                    m.split_face_by_edge(f, v0, v1);
                }
        }
        


        
        
        void console_refit_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("display.refit_trackball");
                return;
            }
            me->refit();
        }

        
        void console_test(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: test");
                return;
            }
            
            me->save_active_mesh();
            me->active_mesh().slit_edges(me->get_vertex_selection());
        }
        
        void console_save(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: save <name->x3d|name->obj> ");
                return;
            }
            string file_name;
            if (args.size()==0) {
                file_name = me->get_file_name();
            }
            else file_name = args[0];
            
            const string extension = file_name.substr(file_name.length()-4,file_name.length());
            if(extension==".obj"){
                obj_save(file_name, me->active_mesh());
                return;
            }
            else if(extension==".off"){
                off_save(file_name, me->active_mesh());
                return;
            }
            else if(extension==".x3d"){
                x3d_save(file_name, me->active_mesh());
                return;
            }
            else if(extension==".bhm") {
                Serialization ser(file_name, std::ios_base::out);
                me->active_mesh().serialize(ser);
                return;
            }
            me->printf("unknown format");
            return;
        }
        
        
        void console_refine_edges(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: refine.split_edges <length>");
                me->printf("splits edges longer than <length>; default is 0.5 times average length");
                return;
            }
            
            me->save_active_mesh();
            
            float thresh = 0.5f;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> thresh;
            }
            
            float avg_length = average_edge_length(me->active_mesh());
            
            refine_edges(me->active_mesh(), thresh * avg_length);
            
            return;
            
        }
        
        void console_cc_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.catmull_clark ");
                me->printf("Does one step of Catmull-Clark subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            cc_split(me->active_mesh(),me->active_mesh());
            cc_smooth(me->active_mesh());
            
            return;
        }
        
        void console_root_cc_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.catmull_clark ");
                me->printf("Does one step of Catmull-Clark subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            rootCC_subdivide(me->active_mesh(),me->active_mesh());
            return;
        }

        
        void console_loop_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.loop");
                me->printf("Does one step of Loop subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            loop_split(me->active_mesh(),me->active_mesh());
            loop_smooth(me->active_mesh());
            
            return;
        }
        
        void console_stitch(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.stitch <rad>");
                me->printf("Stitches faces");
                
                return;
            }
            double r = 1e-30;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> r;
            }
            
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            Vec3d c;
            float rad;
            bsphere(m, c, rad);
            stitch_mesh(me->active_mesh(), r * rad);
            return;
        }
        
        void console_compact(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.compact");
                me->printf("Removes unreferenced vertices");
                
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            m.cleanup();
            return;
        }


        
        void console_remove_val1(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.remove_val1");
                me->printf("Removes valence 1 vertices");
                
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            remove_valence_one_vertices(m);
            return;
        }

        void console_remove_val2(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.remove_val2");
                me->printf("Removes valence 2 vertices");
                
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            remove_valence_two_vertices(m);
            return;
        }

        
        void console_root3_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.root3");
                me->printf("Does one step of sqrt(3) subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            root3_subdivide(me->active_mesh(),me->active_mesh());
            
            return;
        }
        
        
        void console_doosabin_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.doo_sabin ");
                me->printf("Does one step of Doo-Sabin Subdivision");

                return;
            }
            me->save_active_mesh();

            cc_split(me->active_mesh(),me->active_mesh());
            dual(me->active_mesh());

            return;
        }
        
        void console_butterfly_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.butterfly ");
                me->printf("Does one step of Modified Butterfly Subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            butterfly_subdivide(me->active_mesh(),me->active_mesh());
            
            return;
        }
        
        void console_dual(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: dual ");
                me->printf("Produces the dual by converting each face to a vertex placed at the barycenter.");
                return;
            }
            me->save_active_mesh();

            dual(me->active_mesh());

            return;
        }
        
        
        void console_minimize_curvature(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.minimize_curvature <anneal>");
                me->printf("Flip edges to minimize mean curvature.");
                me->printf("If anneal is true, simulated annealing (slow) is used rather than a greedy scheme");
                return;
            }
            me->save_active_mesh();
            
            bool anneal=false;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> anneal;
            }
            
            minimize_curvature(me->active_mesh(), anneal);
            me->post_create_display_list();
            return;
        }
        
        void console_minimize_dihedral(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.minimize_dihedral <iter> <anneal> <use_alpha> <gamma> ");
                me->printf("Flip edges to minimize dihedral angles.");
                me->printf("Iter is the max number of iterations. anneal tells us whether to use ");
                me->printf("simulated annealing and not greedy optimization. use_alpha (default=true) ");
                me->printf("means to use angle and not cosine of anglegamma (default=4) is the power ");
                me->printf("to which we raise the dihedral angle");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1000;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> iter;
            }
            
            bool anneal = false;
            if(args.size() > 1)
            {
                istringstream a0(args[1]);
                a0 >> anneal;
            }
            
            bool use_alpha = true;
            if(args.size() > 2)
            {
                istringstream a0(args[2]);
                a0 >> use_alpha;
            }
            
            float gamma = 4.0f;
            if(args.size() > 3)
            {
                istringstream a0(args[3]);
                a0 >> gamma;
            }
            
            
            minimize_dihedral_angle(me->active_mesh(), iter, anneal, use_alpha, gamma);
            return;
        }
        
        void console_maximize_min_angle(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.maximize_min_angle <thresh> <anneal>");
                me->printf("Flip edges to maximize min angle - to make mesh more Delaunay.");
                me->printf("If the dot product of the normals between adjacent faces < thresh");
                me->printf("no flip will be made. anneal selects simulated annealing rather ");
                me->printf("nthan greedy optimization.");
                return;
            }
            me->save_active_mesh();
            
            float thresh = 0.0f;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> thresh;
            }
            bool anneal = false;
            if(args.size() > 1)
            {
                istringstream a0(args[1]);
                a0 >> anneal;
            }
            maximize_min_angle(me->active_mesh(),thresh,anneal);
            return;
        }
        
        
        void console_optimize_valency(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.valency <anneal> ");
                me->printf("Optimizes valency for triangle meshes. Anneal selects simulated annealing rather than greedy optim.");
                return;
            }
            me->save_active_mesh();
            
            bool anneal = false;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> anneal;
            }
            optimize_valency(me->active_mesh(), anneal);
            return;
        }
        
        
        
        void console_close_holes(MeshEditor* me, const std::vector<std::string> & args)
        {
            int max_size = console_arg(args, 0, 100);

            if(wantshelp(args))
            {
                me->printf("usage: cleanup.close_holes");
                me->printf("This function closes holes. It simply follows the loop of halfvectors which");
                me->printf("enclose the hole and add a face to which they all point.");
                return;
            }
            me->save_active_mesh();
            
            close_holes(me->active_mesh(), max_size);
            return;
        }
        
        void console_reload(MeshEditor* me, const std::vector<std::string> & args)
        {
            string file_name = console_arg(args, 0, me->active_visobj().get_file_name());
            
            if(wantshelp(args))
            {
                me->printf("usage:  load <file>");
                me->printf("(Re)loads the current file if no argument is given, but");
                me->printf("if an argument is given, then that becomes the current file");
                return;
            }
            me->save_active_mesh();

            if(me->reload_active_from_file(file_name))
                me->printf("Loaded %s", file_name.c_str());
            else
                me->printf("failed to load: %s", file_name.c_str());
        }

        
        void console_add_mesh(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  add_mesh <file>");
                me->printf("Loads the file but without clearing the mesh. Thus, the loaded mesh is added to the");
                me->printf("current model.");
                return;
            }
            me->save_active_mesh();
            
            if(!me->add_to_active_from_file(args.size() > 0 ? args[0]:""))
                me->printf("failed to load");
            
            return;
        }
        
        void console_valid(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  validity");
                me->printf("Tests validity of Manifold");
                return;
            }
            VertexSet vs;
            HalfEdgeSet hs;
            FaceSet fs;
            if(find_invalid_entities(me->active_mesh(), vs, hs, fs))
                me->printf("Mesh is valid");
            else {
                me->printf("Mesh is invalid - check console output");
                me->get_vertex_selection() = vs;
                me->get_halfedge_selection() = hs;
                me->get_face_selection() = fs;
            }
            return;
        }
        


        void console_info_all(MeshEditor* me, const std::vector<std::string> & args)
        {
            Vec3d p0_all(FLT_MAX), p7_all(-FLT_MAX);
            for(int i=0;i<me->get_no_meshes();++i)
            {
                Vec3d p0, p7;
                Manifold& m = me->get_mesh(i);
                if (m.no_faces()>0) {
                    bbox(m, p0, p7);
                    p0_all = v_min(p0, p0_all);
                    p7_all = v_max(p7, p7_all);
                }
            }
            
            stringstream bbox_corners;
            bbox_corners << p0_all << " - " << p7_all << endl;
            me->printf("Bounding box corners : %s", bbox_corners.str().c_str());

        }
    
    void console_volume(MeshEditor* me, const std::vector<std::string> & args)
    {
        if(wantshelp(args))
        {
            me->printf("usage:  volume");
            me->printf("Computes mesh volume. Assumes mesh is closed");
            return;
        }
        
        Manifold& m = me->active_mesh();
        for(auto h: m.halfedges())
            if (boundary(m, h)) {
                me->printf("Does not compute: mesh is not watertight");
                return;
            }
        
        Vec3d c;
        float r;
        bsphere(m, c, r);

        double V = 0.0;
        for(auto f: m.faces()) {
            double A = area(m, f);
            Vec3d n = normal(m, f);
            Vec3d p = m.pos(m.walker(f).vertex());
            double h = dot(n, p-c);
            V += A*h;
        }
        V /= 3.0;
        me->printf("Volume = %f", V);
    }
    

        void console_info(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  info");
                me->printf("Provides information about mesh.");
                return;
            }
            Vec3d p0, p7;
            bbox(me->active_mesh(), p0, p7);
            stringstream bbox_corners;
            bbox_corners << p0 << " - " << p7 << endl;
            me->printf("Bounding box corners : %s", bbox_corners.str().c_str());
            map<int,int> val_hist;
            
            Manifold& m = me->active_mesh();
            
            double avg_len = 0;
            for(HalfEdgeID h: m.halfedges()) {
                DebugRenderer::edge_colors[h] = Vec3f(0.3);
                avg_len += length(m,h);
            }
            avg_len /= m.no_halfedges();
            me->printf("Avg. edge length: %f", avg_len);
            for(VertexID v: m.vertices())
            {
                int val = valency(m,v);
                DebugRenderer::vertex_colors[v] = get_color(val);
                ++val_hist[val];
                
                if(val != 4)
                    circulate_vertex_ccw(m, v, static_cast<std::function<void(HalfEdgeID)>>([&](HalfEdgeID h){
                        Walker w = m.walker(h);
                        DebugRenderer::edge_colors[h] = Vec3f(1);
                        DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                        while(valency(m, w.vertex())==4) {
                            w = w.next().opp().next();
                            DebugRenderer::edge_colors[w.halfedge()] = Vec3f(1);
                            DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                        }
                    }));
            }
            map<int, int> ngon_hist;
            for(FaceID f: m.faces()) {
                int ne = no_edges(m, f);
                ++ngon_hist[ne];
                DebugRenderer::face_colors[f] = 0.7*get_color(ne);
            }
            
            me->printf("Valency histogram");
            for(map<int,int>::iterator iter = val_hist.begin(); iter != val_hist.end(); ++iter)
            {
                me->printf("%d, %d", iter->first, iter->second);
            }
            
            me->printf("Ngon histogram");
            for(map<int,int>::iterator iter = ngon_hist.begin(); iter != ngon_hist.end(); ++iter)
            {
                me->printf("%d, %d", iter->first, iter->second);
            }

            
            me->printf("Mesh contains %d faces", me->active_mesh().no_faces());
            me->printf("Mesh contains %d halfedges", me->active_mesh().no_halfedges());
            me->printf("Mesh contains %d vertices", me->active_mesh().no_vertices());
            return;
        }
        
        
        void console_simplify(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: simplify <fraction> <err_thresh> <singular_thresh>");
                me->printf("Performs Garland Heckbert (quadric based) mesh simplification.");
                me->printf("The first argument is the fraction of vertices to keep.");
                me->printf("The second argument is the error threshold relative to bounding sphere radius");
                me->printf("The third argument is the threshold for singular values.");
                return;
            }
            me->save_active_mesh();
            
            float keep_fraction;
            if(args.size() == 0)
            {
                me->printf("you must specify fraction of vertices to keep");
                return;
            }
            istringstream a0(args[0]);
            a0 >> keep_fraction;

            double err_thresh = 1e5;
            if(args.size() >1)
            {
                istringstream a1(args[1]);
                a1 >> err_thresh;
            }

            double singular_thresh = 1e-4;
            if(args.size() >2)
            {
                istringstream a2(args[2]);
                a2 >> singular_thresh;
            }

            cout << "err thresh " << err_thresh << " , sing " << singular_thresh << endl;
        
            
            Vec3d p0, p7;
            Manifold m = me->active_mesh();
            bbox(m, p0, p7);
            Vec3d d = p7-p0;
            float s = d.max_coord();
            Vec3d pcentre = (p7+p0)/2.0;
            for(VertexID v: m.vertices())
            {
                m.pos(v) -= pcentre;
                m.pos(v) /= s;
            }
            cout << "Timing the Garland Heckbert (quadric based) mesh simplication..." << endl;
            Timer timer;
            timer.start();
            
            //simplify
            quadric_simplify(me->active_mesh(),keep_fraction,singular_thresh,err_thresh);
            
            cout << "Simplification complete, process time: " << timer.get_secs() << " seconds" << endl;
            
            for(VertexID v: m.vertices())
            {
                m.pos(v) *= s;
                m.pos(v) += pcentre;
            }
            return;
        }
        
        void console_vertex_noise(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: noise.perturb_vertices <amplitude>");
                me->printf("adds a random vector to each vertex. A random vector in the unit cube is generated and");
                me->printf("to ensure an isotropic distribution, vectors outside the unit ball are discarded.");
                me->printf("The vector is multiplied by the average edge length and then by the amplitude specified.");
                me->printf("If no amplitude is specified, the default (0.5) is used.");
                return;
            }
            me->save_active_mesh();
            
            float avg_length = average_edge_length(me->active_mesh());
            
            float noise_amplitude = 0.5f;
            if(args.size() > 0) {
                istringstream a0(args[0]);
                a0 >> noise_amplitude;
            }
            
            gel_srand(0);
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi){
                Vec3d v;
                do{
                    v = Vec3d(gel_rand(),gel_rand(),gel_rand());
                    v /= (float)(GEL_RAND_MAX);
                    v -= Vec3d(0.5);
                    v *= 2.0;
                }
                while(sqr_length(v) > 1.0);
                
                v *= noise_amplitude;
                v *= avg_length;
                me->active_mesh().pos(*vi) += v;
            }
            return;
        }
        
        void console_perpendicular_vertex_noise(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: noise.perturb_vertices_perpendicular <amplitude>");
                me->printf("adds the normal times a random scalar times amplitude times");
                me->printf("times average edge length to the vertex. (default amplitude=0.5)");
                return;
            }
            me->save_active_mesh();
            
            float avg_length = average_edge_length(me->active_mesh());
            
            float noise_amplitude = 0.5;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> noise_amplitude;
            }
            
            VertexAttributeVector<Vec3d> normals;
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi)
                normals[*vi] = normal(me->active_mesh(), *vi);
            
            gel_srand(0);
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi)
            {
                float rval = 0.5-gel_rand() / float(GEL_RAND_MAX);
                me->active_mesh().pos(*vi) += normals[*vi]*rval*noise_amplitude*avg_length*2.0;
            }
            return;
        }
        
        void console_noisy_flips(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage:  noise.perturb_topology <iter>");
                me->printf("Perform random flips. iter (default=1) is the number of iterations.");
                me->printf("mostly for making nasty synthetic test cases.");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            
            randomize_mesh(me->active_mesh(),  iter);
            return;
        }
        
        void console_laplacian_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  smooth.laplacian <weight> <iter>");
                me->printf("Perform Laplacian smoothing. weight is the scaling factor for the Laplacian.");
                me->printf("default weight = 1.0. Default number of iterations = 1");
                return;
            }
            me->save_active_mesh();
            
            float t=1.0;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> t;
            }
            int iter = 1;
            if(args.size()>1){
                istringstream a0(args[1]);
                a0 >> iter;
            }
            Util::Timer tim;
            tim.start();
            /// Simple laplacian smoothing with an optional weight.
            laplacian_smooth(me->active_mesh(), t, iter);
            cout << "It took "<< tim.get_secs();
            return;
        }
        
        
        void console_taubin_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage:  smooth.taubin <iter>");
                me->printf("Perform Taubin smoothing. iter (default=1) is the number of iterations.");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
            taubin_smooth(me->active_mesh(),  iter);
            
            return;
        }
        
        void console_TAL_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage:  smooth.TAL <weight> <iter>");
                me->printf("Perform Tangential Area weighted Laplacian smoothing. Weight is how strong the");
                me->printf("smoothing should be, and iter (default=1) is the number of iterations.");
                return;
            }
            me->save_active_mesh();
            
            float w=1.0;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> w;
            }
            int iter = 1;
            if(args.size()>1){
                istringstream a0(args[1]);
                a0 >> iter;
            }

            TAL_smoothing(me->active_mesh(), w, iter);
            
            return;
        }

        void console_bilateral_anisotropic_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: smooth.bilateral_anisotropic <iter> <sharpness>");
                me->printf("Smooth normals using a bilateral filter and then update");
                me->printf("geometry. iter (default=1) is the number of iterations. ");
                me->printf("high sharpness (default=4) preserves features better.");
                return;
            }
            me->save_active_mesh();
            
            int iter=1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            double sharpness=4.0;
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> sharpness;
            }

            anisotropic_smooth(me->active_mesh(),  iter, sharpness);
            
            return;
        }
        
        void console_triangulate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  triangulate");
                me->printf("This function triangulates all non triangular faces of the mesh.");
                me->printf("Passing an argument starting with 's' will cause it to split polygons by");
                me->printf("introducing the shortest possible edge. Any other character and it will clip ears.");
                me->printf("You may want to call it after hole closing.");
                return;
            }
            string method = args.empty() ? "S" : args[0];
            for(auto& c:method)
                c = toupper(c);
            
            me->save_active_mesh();
            if(method[0] == 'S') {
                me->printf("Triangulating by shortest edge splits...");
                triangulate(me->active_mesh(), SHORTEST_EDGE);
            }
            else {
                me->printf("Triangulating by ear clipping...");
                triangulate(me->active_mesh(), CLIP_EAR);
            }
            me->active_mesh().cleanup();
            valid(me->active_mesh());
            return;
        }
        
        
        void console_remove_caps(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  cleanup.remove_caps thresh");
                me->printf("Remove caps (triangles with one very big angle). The thresh argument is the fraction of PI to");
                me->printf("use as threshold for big angle. Default is 0.85. Caps are removed by flipping.");
                return;
            }
            me->save_active_mesh();
            
            float t = 0.85f;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> t;
            }
            remove_caps(me->active_mesh(), static_cast<float>(M_PI) *t);
            me->active_mesh().cleanup();
            
            return;
        }
        
        void console_remove_needles(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: cleanup.remove_needles <thresh>");
                me->printf("Removes very short edges by collapse. thresh is multiplied by the median edge length");
                me->printf("to get the length shorter than which we collapse. Default = 0.1");
                return;
            }
            me->save_active_mesh();
            float thresh = console_arg(args, 0, 0.1);
            remove_needles(me->active_mesh(), thresh);
            me->active_mesh().cleanup();
            
            return;
        }
        
        void console_flip_orientation(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: cleanup.flip_orientation");
                me->printf("reorients all faces - flipping normal direction");
                return;
            }
            me->save_active_mesh();
            flip_orientation(me->active_mesh());
            return;
        }

        
        void console_undo(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  undo");
                me->printf("This function undoes one operation. Repeated undo does nothing");
                return;
            }
            me->restore_active_mesh();
            return;
        }
        
        void console_save_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  display.save_trackball");
                me->printf("This function saves the trackball to disk");
                return;
            }
            me->save_ball();
            return;
        }

        void console_load_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  display.load_trackball");
                me->printf("This function loads the trackball from disk");
                return;
            }
            me->load_ball();
            return;
        }

        
    }
    
    
    void MeshEditor::register_console_function(const std::string& name,
                                               const std::function<void(MeshEditor*, const std::vector<std::string>&)>& con_fun,
                                               const std::string& help_txt)
    {
        std::function<void (const std::vector<std::string>&)> f = bind(con_fun, this, placeholders::_1);
        theConsole.reg_cmdN(name, f, help_txt);
    }
    
    void MeshEditor::keyparse(unsigned short key){
        //toggle console with ESC
        if (key == 27)
        {
            console_visible = !console_visible;
            return;
        }
        
        if (console_visible)
        {
            static Timer tim;
            if(key==13)
                tim.start();
            theConsole.keyboard(key);
            if(key == 13)
            {
                active_visobj().post_create_display_list();
                double t = tim.get_secs();
                printf("%f seconds",t);
            }
            return;
        }
        else {

            switch(key) {
                case '\033':
                    console_visible = false;
                    break;
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    active = key - '1';
                    break;
                case ' ':
                    active_visobj().clear_vertex_selection();
                    active_visobj().clear_face_selection();
                    active_visobj().clear_halfedge_selection();
                    post_create_display_list();
                    break;
                case 'f':
                    display_smooth_shading = !display_smooth_shading;
                    post_create_display_list();
                    break;
                case 'w':
                case 'n':
                case 'i':
                case 'r':
                case 't':
                case 'g':
                case 'x':
                case 'a':
                case 'c':
                case 's':
                case 'd':
                case 'C':
                case 'M':
                case 'G':
                    display_render_mode = hotkey_to_string[key];
                    post_create_display_list();
                    break;
                case 'R':
                    theConsole.key_up();
                    theConsole.keyboard(13);
                    post_create_display_list();
                    break;
            }
            
        }
        
    }
    
    void MeshEditor::printf(const char* format, ...)
    {
        //format text
        char buffer[1024];
        va_list args;
        va_start(args, format);
        vsprintf(buffer, format, args);
        va_end(args);
        theConsole.print(buffer);
    }

    void MeshEditor::key_up(){
        if(console_visible)
            theConsole.key_up();
        else
        {
            GLint w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = 0;
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_down(){
        if(console_visible)
            theConsole.key_down();
        else
        {
            GLint w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = NO_MESHES-1;
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_left(){
        if(console_visible)
            theConsole.key_left();
        else
        {
            GLint w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = max(0, active-1);
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_right(){
        if(console_visible)
            theConsole.key_right();
        else
        {
            GLint w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = min(NO_MESHES-1, active+1);
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_home(){
        theConsole.key_home();
    }
    void MeshEditor::key_end(){
        theConsole.key_end();
    }
   
    bool MeshEditor::listen_commands() {
        static Timer tim;
        tim.start();
        if(theConsole.listen_commands()) {
            double t = tim.get_secs();
            printf("%f seconds",t);
// Takes too long: post_create_display_list();
            return true;
        }
        return false;
    }
    
    void MeshEditor::grab_ball(TrackBallAction action, const CGLA::Vec2i& pos){
        active_view_control().grab_ball(action, pos);
    }
    void MeshEditor::roll_ball(const CGLA::Vec2i& pos){
        active_view_control().roll_ball(pos);
    }
    void MeshEditor::release_ball(){
        active_view_control().release_ball();
    }
    bool MeshEditor::try_spinning_ball(){
        return active_view_control().try_spin();
    }
    
    void MeshEditor::save_ball() {
        ofstream ofs("trackball.bin", ios_base::binary);
        active_view_control().save(ofs);
        
    }
    void MeshEditor::load_ball() {
        ifstream ifs("trackball.bin", ios_base::binary);
        active_view_control().load(ifs);
    
    }

    
    bool MeshEditor::grab_mesh(const CGLA::Vec2i& pos)
    {
        if(depth_pick(pos[0], pos[1], depth))
        {
            dragging = true;
            mouse_x = pos[0];
            mouse_y = pos[1];
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Manifold& m = active_mesh();
            active_visobj().save_old();
            Vec3d c;
            float r;
            bsphere(m, c, r);
            for(auto vid : m.vertices())
            {
                double l = sqr_length(p0-m.pos(vid));
                weight_vector[vid] = exp(-l/(brush_size*r*r));
            }
            return true;
        }
        return false;
    }
    
    bool MeshEditor::drag_mesh(const CGLA::Vec2i& pos)
    {
        auto deform_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, Vec3d& v)
        {
            const Manifold& m_old = active_visobj().mesh_old();
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m_old.vertices())
                new_pos[vid] = m_old.pos(vid) + weight_vector[vid] * v;
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto inflate_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * (0.25 * laplacian(m, vid) + (r*0.0025)*normal(m, vid));
            }
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto smooth_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * 0.5 * laplacian(m, vid);
            }
            m.positions_attribute_vector() = new_pos;
        };

        auto paint_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            auto& col_map = active_visobj().get_color_field_attrib_vector();
            
            /// This is inelegant, but we need to know if the damn thing is initialized.
            if(std::isnan(col_map[*m.vertices().begin()][0])) {
                cout << "col_map.size " << col_map.size() << endl;
                for(auto vid: m.vertices())
                    col_map[vid] = Vec3d(0);
            }
            
            double support_radius = r * brush_size;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = length(p0-p);
                if(l<support_radius) {
                    double t = l/support_radius;
                    double wgt = 1.0 - (3*t*t - 2 * t*t*t);
                    col_map[vid] += 0.1*wgt*Vec3d(paint_color);
                }
            }
        };

        
        if(dragging)
        {
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Vec3d p1 = screen2world(pos[0], pos[1], depth);
            Vec3d v = p1-p0;
            Manifold& m = active_mesh();
            auto vset = active_visobj().get_vertex_selection();
            auto hset = active_visobj().get_halfedge_selection();
            auto fset = active_visobj().get_face_selection();
            
            if(!vset.empty())
                for(auto vid : vset)
                    m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
            else if(!hset.empty())
                for(auto h : hset) {
                    Walker w = m.walker(h);
                    auto vid0 = w.vertex();
                    auto vid1 = w.opp().vertex();
                    m.pos(vid0) = active_visobj().mesh_old().pos(vid0) + v;
                    m.pos(vid1) = active_visobj().mesh_old().pos(vid1) + v;
                }
            else if(!fset.empty())
                for(auto f: fset)
                {
					circulate_face_ccw(m, f, std::function<void(VertexID)>( [&](VertexID vid){
                        m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
                    }) );
                }
            else {
                if(string(brush_type) == "smooth")
                    smooth_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "inflate")
                    inflate_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "deform")
                    deform_mesh(m, weight_vector, v);
                else if(string(brush_type) == "paint") {
                    float d;
                    if(depth_pick(pos[0], pos[1], d))
                        paint_mesh(m, weight_vector, screen2world(pos[0], pos[1], d));
                }
            }
            
            post_create_display_list();
            return true;
        }
        return false;
    }
    
    void MeshEditor::release_mesh()
    {
        dragging = false;
    }
    
    
    
    
//    mutex parallel_work;
    void MeshEditor::display(int scale){

        // Clear screen.
        glClearColor(1, 1, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        //        parallel_work.lock();
        //        active_visobj().post_create_display_list();
        
        
        // Display object.
        active_visobj().display(display_render_mode, theConsole, display_smooth_shading, display_gamma);
        
        // Draw console
        if(console_visible)
        {
            glUseProgram(0);
            theConsole.display(scale);
        }
        
        // Static variable controlling whether we render shadow at all.
        static Console::variable<int> shadow_enable(0);
        shadow_enable.reg(theConsole, "display.shadow.enable", "");
        if(shadow_enable) {
            // Static variables that control the shadow display created and linked to console
            static Console::variable<float> zenith(1.571);
            static Console::variable<float> azimuth(0);
            static Console::variable<float> shadow_alpha(0.3);
            azimuth.reg(theConsole, "display.shadow.azimuth", "");
            zenith.reg(theConsole, "display.shadow.zenith", "");
            shadow_alpha.reg(theConsole, "display.shadow.shadow_alpha", "");

            
            // Shadow buffer (really a frame buffer object)
            static ShadowBuffer sb;
            
            // String containing fragment program for shadow rendering
            static string shadow_shdr_fp = "#version 120\n"
            "    uniform sampler2DShadow shadow_map;\n"
            "    uniform mat4 Mat;\n"
            "   uniform float shadow_alpha;\n"
            "    varying vec4 ep;\n"
            "    void main()\n"
            "    {\n"
            "        vec4 light_pos =  Mat * ep;\n"
            "        light_pos.z = max(0.001,min(0.999,light_pos.z-0.003));\n"
            "        if(light_pos.x <=0 || light_pos.x >=1 || light_pos.y <=0 || light_pos.y>=1) gl_FragColor = vec4(1);"
            "   else      gl_FragColor= vec4(1.0-shadow_alpha)+shadow_alpha*0.25*(shadow2D(shadow_map, light_pos.xyz+vec3(0.001,-0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz-vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(-0.001,0.001,0)));\n"
            "    }\n";
            
            // Shader program for shadow rendering is compiled and linked.
            static GLuint prog = 0;
            if(!prog)
            {
                GLuint vp = create_glsl_shader(GL_VERTEX_SHADER,"#version 120\n varying vec4 ep; void main(){ep = gl_Vertex; gl_Position = ftransform();}");
                GLuint fp = create_glsl_shader(GL_FRAGMENT_SHADER, shadow_shdr_fp);
                prog = glCreateProgram();
                
                glAttachShader(prog, vp);
                glAttachShader(prog, fp);
                glLinkProgram(prog);
            }
            
            
            // Setup OpenGL state for lighting - used when rendering 3D object casting shadow
            Vec4f lpos = Vec4f(cos(azimuth)*cos(zenith),sin(zenith),sin(azimuth)*cos(zenith),0);
            glLightfv(GL_LIGHT0, GL_POSITION, lpos.get());
            Vec4f mamb(.8,.8,.8,0);
            Vec4f lamb(.4,.4,.5,0);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mamb.get());
            glLightfv(GL_LIGHT0, GL_AMBIENT, lamb.get());

            // Get old viewport dimensions. Setup rendering to FBO and set viewport
            // dimensions for rendering shadows.
            GLint viewp[4];
            glGetIntegerv(GL_VIEWPORT, viewp);
            sb.enable();
            glViewport(0, 0, 1024, 1024);
            
            // Setup object transformations using old school GL.
            Mat4x4f m;
            float r = active_visobj().get_bsphere_radius();
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            glOrtho(-2*r, 2*r, -2*r, 2*r, -2*r, 2*r);
            gluLookAt(0.1*r*cos(azimuth)*cos(zenith),0.1*r*sin(zenith),0.1*r*sin(azimuth)*cos(zenith), 0,0,0, 0,1,0);
            
            // Copy the transformation matrix to user code.
            glGetFloatv(GL_PROJECTION_MATRIX, m.get());
            
            // Draw the object in light space.
            draw(active_visobj().mesh());
            
            // Restore transformation matrices.
            glPopMatrix();
            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();
            
            // Restore the usual on screen framebuffer.
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
            glViewport(0, 0, viewp[2], viewp[3]);
            
            // Use shadow: bind the shadow texture to unit 0.
            sb.bind_textures(0);
            
            // Use the shadow rendering shader program.
            glUseProgram(prog);
            active_visobj().view_control().set_gl_modelview();
            glUniform1i(glGetUniformLocation(prog, "shadow_map"), 0);
            glUniform1f(glGetUniformLocation(prog, "shadow_alpha"), shadow_alpha);
            m = translation_Mat4x4f(Vec3f(0.5)) * scaling_Mat4x4f(Vec3f(0.5)) * transpose(m);
            glUniformMatrix4fv(glGetUniformLocation(prog, "Mat"), 1, 1, m.get());
            
            
            // Setup blending such that the shadow is alpha blended with model.
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ZERO, GL_SRC_COLOR);

            // Draw ground plane for shadows.
            Vec3d p0, p7;
            bbox(active_visobj().mesh(), p0, p7);
            glBegin(GL_QUADS);
            glVertex3f(-100*r, p0[1],-100*r);
            glVertex3f(-100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],-100*r);
            glEnd();
            
            // Draw model again ... just to add shadow.
            draw(active_visobj().mesh());
            
            // Disable blending and shader program.
            glDisable(GL_BLEND);
            glUseProgram(0);
        }
        //        parallel_work.unlock();
        //        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    void MeshEditor::reshape(int w, int h) {
        for(VisObj& v : vo)
            v.view_control().reshape(w, h);
    }
    
    Vec2i MeshEditor::shape() {
        return vo[active].view_control().shape();
    }

    void MeshEditor::init() {
        glewInit();
        
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT, vp);
        for(VisObj& vis_obj : vo)
            vis_obj.view_control().reshape(vp[2], vp[3]);
        
//        glEnable(GL_CULL_FACE);
//        glCullFace(GL_BACK);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glClearColor(1,1,0, 0.f);
        glColor4f(1.0f, 1.0f, 1.0f, 0.f);
        float material[4] = {1,1,1,1};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material);
        glEnable(GL_DEPTH_TEST);
        
        register_console_function("quit", console_quit,"");
        register_console_function("exit", console_quit,"");
        register_console_function("bye", console_quit,"");
        
        register_console_function("merge_with", console_merge, "");
        register_console_function("clear_mesh", console_clear, "");

        register_console_function("help.controls", console_list_controls, "");
        register_console_function("simplify", console_simplify,"");
        register_console_function("ridge_lines", console_ridge_lines,"");
        
        register_console_function("smooth.laplacian", console_laplacian_smooth,"");
        register_console_function("smooth.taubin", console_taubin_smooth,"");
        register_console_function("smooth.TAL", console_TAL_smooth,"");
        register_console_function("smooth.bilateral_anisotropic", console_bilateral_anisotropic_smooth ,"");
        
        register_console_function("optimize.valency", console_optimize_valency,"");
        register_console_function("optimize.minimize_dihedral_angles", console_minimize_dihedral,"");
        register_console_function("optimize.minimize_curvature", console_minimize_curvature,"");
        register_console_function("optimize.maximize_min_angle", console_maximize_min_angle,"");
        register_console_function("load_mesh", console_reload,"");
        register_console_function("add_mesh", console_add_mesh,"");
        
        register_console_function("cleanup.stitch", console_stitch,"");
        register_console_function("cleanup.remove_val1", console_remove_val1, "");
        register_console_function("cleanup.remove_val2", console_remove_val2, "");
        register_console_function("cleanup.flip_orientation", console_flip_orientation,"");
        register_console_function("cleanup.remove_caps", console_remove_caps,"");
        register_console_function("cleanup.remove_needles", console_remove_needles,"");
        register_console_function("cleanup.close_holes", console_close_holes,"");
        register_console_function("cleanup.compact", console_compact,"");

        register_console_function("triangulate", console_triangulate,"");
        register_console_function("dual", console_dual,"");
        register_console_function("refine.split_edges", console_refine_edges,"");
        
        register_console_function("subdivide.catmull_clark", console_cc_subdivide,"");
        register_console_function("subdivide.rootcc", console_root_cc_subdivide,"");
        register_console_function("subdivide.loop", console_loop_subdivide,"");
        register_console_function("subdivide.root3", console_root3_subdivide,"");
        register_console_function("subdivide.doo_sabin", console_doosabin_subdivide,"");
        register_console_function("subdivide.butterfly", console_butterfly_subdivide,"");
        register_console_function("save_mesh", console_save,"");
        register_console_function("noise.perturb_vertices", console_vertex_noise,"");
        register_console_function("noise.perturb_vertices_perpendicular", console_perpendicular_vertex_noise,"");
        register_console_function("noise.perturb_topology", console_noisy_flips,"");
        
        
        register_console_function("align_with", console_align,"");
        register_console_function("undo", console_undo,"");
        
        register_console_function("validity", console_valid,"");
        register_console_function("info", console_info,"");
        register_console_function("volume", console_volume,"");
        register_console_function("info.all_meshes", console_info_all,"");
        
                
        register_console_function("display.refit_trackball", console_refit_trackball, "Resets trackball");
        register_console_function("display.save_trackball", console_save_trackball, "Saves trackball to disk");
        register_console_function("display.load_trackball", console_load_trackball, "Load trackball to disk");
        
        register_console_function("transform.scale", console_scale, "Scale mesh");
        register_console_function("transform.rotate", console_rotate, "Rotate mesh");
        register_console_function("transform.translate", console_translate, "Translate mesh");
        
        register_console_function("edit.selected.merge_1_ring", console_merge_1_ring, "Merge 1-ring of selected vertices");
        
        register_console_function("edit.selected.split_edge", console_split_edge , "");
        register_console_function("edit.selected.collapse_edge", console_collapse_edge, "");
        register_console_function("edit.selected.flip_edge", console_flip_edge, "");
        register_console_function("edit.selected.dissolve_edge", console_dissolve_edge, "");
        register_console_function("edit.selected.stellate_face", console_stellate_face, "");
        register_console_function("edit.selected.split_face", console_split_face, "");
        register_console_function("edit.selected.delete", console_delete, "");
        register_console_function("edit.selected.triangulate_face", console_triangulate_face, "");
        register_console_function("edit.selected.bridge_faces", console_bridge_faces, "");

        register_console_function("test", console_test, "Test some shit");
        
        selection_mode.reg(theConsole, "selection.mode", "The selection mode. 0 = vertex, 1 = halfedge, 2 = face");
        active.reg(theConsole, "active_mesh", "The active mesh");
        display_render_mode.reg(theConsole, "display.render_mode", "Display render mode");
        brush_size.reg(theConsole, "brush.size", "Size of brush used for editing");
        brush_type.reg(theConsole, "brush.type", "Smooth, deform, inflate, or paint");
        paint_color.reg(theConsole, "brush.color", "The color of the brush used in paint mode");
        
        display_smooth_shading.reg(theConsole, "display.smooth_shading", "1 for smooth shading 0 for flat");
        display_gamma.reg(theConsole, "display.gamma", "The gamma setting for the display");

        theConsole.print("Welcome to MeshEdit");
        theConsole.newline();
    }
    
    bool MeshEditor::add_file(const std::string& str)
    {   
        while (active_mesh().no_vertices()>0 && active<NO_MESHES)
            active  = active + 1;
        if(active == NO_MESHES)
            active = 0;
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }

    bool MeshEditor::reload_active_from_file(const std::string& str)
    {
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }
    
    bool MeshEditor::add_to_active_from_file(const std::string& str)
    {
        if(active_visobj().add_mesh(str)) {
            active_visobj().post_create_display_list();
            return  true;
        }
        return false;
    }

    
}
