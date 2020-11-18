//
//  MeshEditor.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/10/13.
//
//

#ifndef __GEL__MeshEditor__
#define __GEL__MeshEditor__

#include <sstream>
#include <mutex>
#include <string>
#include <GEL/GLGraphics/Console.h>
#include <GEL/GLGraphics/VisObj.h>
#include <GEL/GLGraphics/GLViewController.h>


namespace GLGraphics {
    extern std::mutex parallel_work;
    
    
    template<typename T>
    T console_arg(const std::vector<std::string> & args, int num, T dflt)
    {
        if(args.size() > num){
            std::istringstream iss(args[num]);
            T a;
            iss >> a;
            return a;
        }
        return dflt;
    }

    
    class MeshEditor
    {
        bool console_visible = false;
        bool dragging = false;
        int mouse_x, mouse_y;
        float depth;
        HMesh::VertexAttributeVector<float> weight_vector;
        HMesh::VertexAttributeVector<CGLA::Vec3d> orig_pos;

        static const int NO_MESHES = 25;
        std::array<VisObj,NO_MESHES> vo;
        
    public:
        VisObj& active_visobj() { return vo[active]; }
        const VisObj& active_visobj() const { return vo[active]; }
    private:
        GLViewController& active_view_control() {
            return active_visobj().view_control();
        }

        Console theConsole;
        Console::variable<std::string> brush_type;
        Console::variable<int> selection_mode;
        Console::variable<int> active;
        Console::variable<std::string> display_render_mode;
        Console::variable<float> brush_size;
        Console::variable<CGLA::Vec3d> paint_color;
        Console::variable<int> display_smooth_shading;
        Console::variable<float> display_gamma;

    public:
        MeshEditor():
        brush_type("smooth"),
        selection_mode(0),
        active(0),
        display_render_mode("iso"),
        brush_size(0.01),
        display_smooth_shading(true),
        display_gamma(2.2) {}

        /// Initialize the mesh editor. Do this only when OpenGL state is available.
        void init();
        
        bool select(const CGLA::Vec2i& pos) {
            switch(selection_mode){
                case 0: return active_visobj().select_vertex(pos);
                case 1: return active_visobj().select_halfedge(pos);
                case 2: return active_visobj().select_face(pos);
            }
            return false;
        }
        
         HMesh::VertexSet& get_vertex_selection()  {
            return active_visobj().get_vertex_selection();
        }
        
        HMesh::FaceSet& get_face_selection()  {
            return active_visobj().get_face_selection();
        }

        HMesh::HalfEdgeSet& get_halfedge_selection()  {
            return active_visobj().get_halfedge_selection();
        }

        /// Returns the number of mesh objects
        int get_no_meshes() const {return NO_MESHES;}

        
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
        CGLA::Vec2i shape();
        
        
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
        void refit(const CGLA::Vec3d& c, double r) {
            active_visobj().refit(c,r);
        }
        void refit() {
            active_visobj().refit();
        }

        /// Returns the name of the file whence the active mesh was loaded.
        std::string& get_file_name() {return active_visobj().get_file_name();}

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
        
        bool listen_commands();

        /// Returns a reference to active mesh.
        HMesh::Manifold& active_mesh() { return active_visobj().mesh(); }
        
        
        void set_active(int i) {
            active = std::min((NO_MESHES-1), (std::max(0, i)));
        }
        
        int get_active_no() const {
            return active;
        }

        /// Returns a reference to mesh i
        HMesh::Manifold& get_mesh(int i) { return vo[i].mesh(); }
        
        /// Returns a reference to visobj i
        GLGraphics::VisObj& get_visobj(int i) { return vo[i];}

       /// Add a file to the next empty slot.
        bool add_file(const std::string& str);
        
        /// Load the mesh given as argument to the current slot.
        bool reload_active_from_file(const std::string& str);
        
        /// Load the mesh but without clearing, effectively combining it with existing mesh.
        bool add_to_active_from_file(const std::string& str);
        
    };
}

#endif /* defined(__GEL__MeshEditor__) */
