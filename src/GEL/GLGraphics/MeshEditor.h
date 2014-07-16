//
//  MeshEditor.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/10/13.
//
//

#ifndef __GEL__MeshEditor__
#define __GEL__MeshEditor__

#include <mutex>
#include <string>
#include "../GLGraphics/Console.h"
#include "../GLGraphics/VisObj.h"
#include "../GLGraphics/GLViewController.h"


namespace GLGraphics {
    extern std::mutex parallel_work;
    
    class MeshEditor
    {
        bool console_visible = false;
        
        bool dragging = false;
        int mouse_x, mouse_y;
        float depth;
        HMesh::VertexAttributeVector<float> weight_vector;
        HMesh::VertexAttributeVector<CGLA::Vec3d> orig_pos;

        static const int NO_MESHES = 9;
        std::array<VisObj,NO_MESHES> vo;
        
    public:
        VisObj& active_visobj() { return vo[active]; }
        const VisObj& active_visobj() const { return vo[active]; }
    private:
        GLViewController& active_view_control() {
            return active_visobj().view_control();
        }

        Console theConsole;
        Console::variable<int> active;
        Console::variable<std::string> display_render_mode;
        Console::variable<float> brush_size;
        Console::variable<int> display_smooth_shading;
        Console::variable<float> display_gamma;

    public:
        MeshEditor():active(0), display_render_mode("normal"), brush_size(0.01), display_smooth_shading(true),
        display_gamma(2.2) {}

        /// Initialize the mesh editor. Do this only when OpenGL state is available.
        void init();
        
        bool select_vertex(const CGLA::Vec2i& pos) {
            return active_visobj().select_vertex(pos);
        }
        
         HMesh::VertexAttributeVector<int>& get_vertex_selection()  {
            return active_visobj().get_vertex_selection();
        }
        

        
        /** Tests whether the position passed as argument is on the mesh (return true) 
         or background (false). This function also retains the 3D unprojection of the 
         grabbed position. */
        bool grab_mesh(const CGLA::Vec2i& pos);
        
        /** Provided grab_mesh has been called and returned true, drag_mesh computes a 
         vector to the new position given as argument and moves a small subset of the mesh
         according to the vector from grabbed to dragged position. */
        bool drag_mesh(const CGLA::Vec2i& pos);
        
        /** Releases the mesh. We are no longer dragging. */
        void release_mesh();

        // GLViewController stuff
        void reshape(int w, int h);
        void grab_ball(TrackBallAction action, const CGLA::Vec2i& pos);
        void roll_ball(const CGLA::Vec2i& pos);
        void release_ball();
        bool try_spinning_ball();
        void save_ball();
        void load_ball();
        
        /// Align means that we sync view controllers.
        void align(int src, int dst)
        {
            vo[dst].view_control() =
            vo[src].view_control();
        }
        
        /// Make sure the object fits in the window.
        void refit() {
            active_visobj().refit();
        }
        
        /// Returns the name of the file whence the active mesh was loaded.
        const std::string& file_name() const {return active_visobj().file_name();}

        // Get mesh and mesh state.
        void save_active_mesh() {active_visobj().save_old();}
        void restore_active_mesh() {active_visobj().restore_old();}

        // Display functions ------------
        
        /// Notify the visualization object that we need to regenerate the display list.
        void post_create_display_list() {active_visobj().post_create_display_list();}
        
        /// Render the active mesh.
        void display(int scale = 1);

        
        // Console functions --------------
        
        void register_console_function(const std::string& name,
                                       const std::function<void(MeshEditor*, const std::vector<std::string>&)>& con_fun,
                                       const std::string& help_txt);

        void printf(const char* format, ...);
        void keyparse(unsigned short c);
        void key_up();
        void key_down();
        void key_left();
        void key_right();
        void key_home();
        void key_end();

        /// Returns a reference to active mesh.
        HMesh::Manifold& active_mesh() { return active_visobj().mesh(); }

        /// Returns a reference to mesh i
        HMesh::Manifold& get_mesh(int i) { return vo[i].mesh(); }

       /// Add a file to the next empty slot.
        bool add_file(const std::string& str);
        
        /// Load the mesh given as argument to the current slot.
        bool reload_active_from_file(const std::string& str);
        
        /// Load the mesh but without clearing, effectively combining it with existing mesh.
        bool add_to_active_from_file(const std::string& str);
        
        void harmonics_analyze_mesh() {active_visobj().harmonics_analyze();}
        void harmonics_reset_shape() {active_visobj().harmonics_reset_shape();}
        void harmonics_partial_reconstruct(int E0, int E1, float scale) {active_visobj().harmonics_partial_reconstruct(E0, E1, scale);}

    };
}

#endif /* defined(__GEL__MeshEditor__) */
