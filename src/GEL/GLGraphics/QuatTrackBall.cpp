/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <iostream>
#include "../GL/glew.h"
#include "../CGLA/CGLA.h"
#include "QuatTrackBall.h"
#include "GLViewController.h"


using namespace std;
using namespace CGLA;

namespace GLGraphics
{
    
    QuatTrackBall::QuatTrackBall(const Vec3f& _centre,
                                 float _eye_dist,
                                 unsigned _width,
                                 unsigned _height):
	centre(_centre), width(_width), height(_height),
	eye_dist(_eye_dist)
	{
        // This size should really be based on the distance from the center of
        // rotation to the point on the object underneath the mouse.  That
        // point would then track the mouse as closely as possible.  This is a
        // simple example, though, so that is left as an exercise.
        ballsize = 2.0f;
        screen_centre = Vec2i(width/2, height/2);
        qrot = Quatf(0.0, 0.0, 0.0, 1.0);
        qinc = Quatf(0.0, 0.0, 0.0, 1.0);
        trans = Vec2f(0.0, 0.0);
    }
    
    void QuatTrackBall::grab_ball(TrackBallAction act, const Vec2i& v)
    {
        if(v[0] < 0 || v[0] >= static_cast<int>(width)
           || v[1] < 0 || v[1] >= static_cast<int>(height))
            return;
        
        set_position(scalePoint(v));
        current_action = act;
    }
    
    void QuatTrackBall::roll_ball(const Vec2i& v)
    {
        if(v[0] < 0 || v[0] >= static_cast<int>(width)
           || v[1] < 0 || v[1] >= static_cast<int>(height))
            return;
        
        Vec2f w = scalePoint(v);
        
        switch (current_action)
        {
			case ROTATE_ACTION:
                rotate(w);
                break;
                
			case PAN_ACTION:
                pan(w);
                break;
                
			case ZOOM_ACTION:
                zoom(w);
                break;
            case NO_ACTION:
            default:
                break;
        }
        last_pos = w;
    }
    
    // Call this when the user does a mouse down.
    // Stop the trackball glide, then remember the mouse
    // down point (for a future rotate, pan or zoom).
    void QuatTrackBall::set_position(const Vec2f& _last_pos)
    {
        stop_spin();
        last_pos = _last_pos;
    }
    
    // Rotationaly spin the trackball by the current increment.
    // Use this to implement rotational glide.
    void QuatTrackBall::do_spin()
    {
        qrot = qrot*qinc;
    }
    
    // Cease any rotational glide by zeroing the increment.
    void QuatTrackBall::stop_spin()
    {
        qinc.set(0.0, 0.0, 0.0, 1.0);
    }
    
    void QuatTrackBall::rotate(const Vec2f& new_v)
    {
        calcRotation(new_v);
        do_spin();
    }
    
    void QuatTrackBall::pan(const Vec2f& new_v)
    {
        trans += (new_v - last_pos) * Vec2f(eye_dist);
    }
    
    void QuatTrackBall::zoom(const Vec2f& new_v)
    {
        eye_dist += (new_v[1] - last_pos[1]) * eye_dist;
        eye_dist = max(eye_dist, 0.01f);
    }
    
    void QuatTrackBall::calcRotation(const Vec2f& new_pos)
    {
        // Check for zero rotation
        if (new_pos == last_pos)
            qinc = Quatf(0.0f, 0.0f, 0.0f, 1.0f);
        else
        {
            // Form two vectors based on input points, find rotation axis
            Vec3f p1 = Vec3f(new_pos[0], new_pos[1], projectToSphere(new_pos));
            Vec3f p2 = Vec3f(last_pos[0], last_pos[1], projectToSphere(last_pos));
            qinc.make_rot(normalize(p1), normalize(p2));
            /*
             Vec3f q = cross(p1, p2);		// axis of rotation from p1 and p2
             float L = sqrt(1.0f-dot(q,q) / (dot(p1,p1) * dot(p2,p2)));
             
             q.normalize();				// q' = axis of rotation
             q *= sqrt((1 - L)/2);	// q' = q' * sin(phi)
             
             qinc.set(q[0],q[1],q[2],sqrt((1 + L)/2));
             */
        }
    }
    
    // Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
    // if we are away from the center of the sphere.
    float QuatTrackBall::projectToSphere(const Vec2f& v)
    {
#ifndef M_SQRT_2
        const double M_SQRT_2 = 0.707106781187;
#endif
        
        float d = v.length();
        float t = ballsize*M_SQRT_2;
        float z;
        
        // Inside sphere
        if(d < ballsize)
            z = sqrt(ballsize*ballsize - d*d);
        else if(d < t)
            z = 0.0;
        // On hyperbola
        else
            z = t*t/d;
        
        return z;
    }
    
    // Scales integer point to the range [-1, 1]
    Vec2f QuatTrackBall::scalePoint(const Vec2i& v) const
    {
        Vec2f w(v[0], v[1]);
        w -= Vec2f(screen_centre);
        w /= Vec2f(width,height);
        w = CGLA::v_min(Vec2f(1.0f), CGLA::v_max(Vec2f(-1), 2*w));
        return w;
    }
    
    void QuatTrackBall::get_view_param(Vec3f& eye, Vec3f& _centre, Vec3f& up) const
    {
        up  = qrot.apply_unit(Vec3f(0,1,0));
        Vec3f right = qrot.apply_unit(Vec3f(1,0,0));
        _centre = centre - up * trans[1] - right * trans[0];
        eye = qrot.apply_unit(Vec3f(0,0,1)*eye_dist) + _centre;
    }
    
    
    // Modify the current gl matrix by the trackball rotation and translation.
    void QuatTrackBall::set_gl_modelview() const
    {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        Vec3f eye;
        Vec3f _centre;
        Vec3f up;
        get_view_param(eye, _centre, up);
        gluLookAt(eye[0], eye[1], eye[2],
                  _centre[0], _centre[1], _centre[2], 
                  up[0],up[1],up[2]);
    }
    
    bool QuatTrackBall::is_spinning() const
    {
        static const Quatf null_quat(0,0,0,1);
        if(!(qinc == null_quat))
            return true;
        return false;
    }
}
