/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../GL/glew.h"
#include "GLViewController.h"
#include "../CGLA/Mat3x3f.h"

using namespace std;
using namespace CGLA;

namespace GLGraphics
{
    
    void GLViewController::reset_projection()
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(FOV_DEG, aspect, znear, zfar);
        glMatrixMode(GL_MODELVIEW);
    }
    
    GLViewController::GLViewController()
    : FOV_DEG(53),
    WINX(500), WINY(500),
    aspect(500/float(500)),
    button_down(false),
    spin(false),
    ball(CGLA::Vec3f(0.0,0.0,0.0), 0.0, 500, 500)
    {
        znear = 0.01f*0.0;
        zfar  = 3*0.0;
        //    reset_projection();
    }
    
    
    GLViewController::GLViewController(int _WINX, int _WINY, const CGLA::Vec3f& centre, float rad)
    : FOV_DEG(53),
    WINX(_WINX), WINY(_WINY),
    aspect(WINX/float(WINY)),
    button_down(false),
    spin(false),
    ball(centre, rad, WINX, WINY)
    {
        znear = 0.01f*rad;
        zfar  = 3*rad;
        // reset_projection();
    }
    
    void GLViewController::grab_ball(TrackBallAction action, const CGLA::Vec2i& pos)
    {
        ball.grab_ball(action,pos);
        if(action==ZOOM_ACTION)
            set_near_and_far();
        
        spin = false;
        button_down = true;
        last_action = action;
    }
    
    void GLViewController::roll_ball(const CGLA::Vec2i& pos)
    {
        static Vec2i old_pos = pos;
        Vec2f dir = Vec2f(pos-old_pos);
        float len = dir.length();
        if (len < TINY)
            return;
        
        ball.roll_ball(pos);
        if(last_action==ZOOM_ACTION)
            set_near_and_far();
        
        spin = len>=1.1f;
        old_pos = pos;
    }
    
    
    void GLViewController::release_ball()
    {
        ball.release_ball();
        if(last_action==ZOOM_ACTION)
            set_near_and_far();
    }
    
    bool GLViewController::try_spin()
    {
        if(spin && !ball.is_grabbed())
        {
            ball.do_spin();
            return true;
        }
        return false;
    }
    
    void GLViewController::set_gl_modelview()
    {
        reset_projection();
        ball.set_gl_modelview();
    }
    
    
    void GLViewController::reshape(int W, int H)
    {
        WINX = W;
        WINY = H;
        aspect = WINX/static_cast<float>(WINY);
        glViewport(0,0,WINX,WINY);
        reset_projection();
        ball.set_screen_window(WINX, WINY);
    }
    
    void GLViewController::set_near_and_far()
    {
        float rad = ball.get_eye_dist();
        znear = 0.01f*rad;
        zfar = 3*rad;
        reset_projection();
    }
    
    void GLViewController::set_view_param(const Vec3f& e, const Vec3f& c, const Vec3f& u)
    {
        // native viewing direction is the negative z-axis
        // while right is the x-axis and up is the y-axis
        Vec3f view = c - e;
        float eye_dist = length(view);
        view /= eye_dist;
        Vec3f right = normalize(cross(view, u));
        Vec3f up = cross(right, view);
        Mat3x3f rot(right, up, -view);
        rot = transpose(rot); // since matrix is row-major
        
        // convert the change-of-basis matrix to a quaternion
        Quatf qrot;
        qrot.make_rot(rot);
        set_rotation(qrot);
        set_centre(c);
        set_eye_dist(eye_dist);
    }
    
    bool GLViewController::load(std::ifstream& ifs)
    {
        if(ifs)
        {
            ifs.read(reinterpret_cast<char*>(this), sizeof(GLViewController));
            reset_projection();
            ball.set_screen_window(WINX, WINY);
            return true;
        }
        return false;
    }
    
    bool GLViewController::save(std::ofstream& ofs) const
    {
        if(ofs)
        {
            ofs.write(reinterpret_cast<const char*>(this), sizeof(GLViewController));
            return true;
        }
        return false;
    }
}
