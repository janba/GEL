#ifndef CAMERA_H
#define CAMERA_H

#include <GEL/GL/glew.h>

#include "GEL/CGLA/Vec2f.h"
#include "GEL/CGLA/Vec3f.h"

#include "GEL/Geometry/Ray.h"

const float NEAR_PLANE = 1.0e-2f;
const float FAR_PLANE = 1.0e8f;

class Camera
{
public:
  
  Camera(const CGLA::Vec3f& _eye,     // Eye point (camera position)
	 const CGLA::Vec3f& _focus,   // Focus point
	 const CGLA::Vec3f& _up,      // Up vector
	 float fd)                    // Focal distance
    : eye(_eye), focus(_focus), up(_up), focal_dist(fd)
  {
#ifndef M_1_PI
    const double M_1_PI = 0.318309886184;
#endif

	line_of_sight = focus - eye;

    // Calculate view plane normal and basis
    vp_normal = normalize(line_of_sight);
    assert(dot(vp_normal, normalize(up)) < 0.99999999);
    vp_axes[0] = normalize(cross(vp_normal, up));
    vp_axes[1] = normalize(cross(vp_axes[0], vp_normal));

    // Account for focal distance
    vp_normal *= focal_dist;
  
    // Calculate field of view (using the pinhole camera model)
    float tmp = atan(1.0/(2.0 * focal_dist));
    fov = 360.0 * M_1_PI * tmp;
  }

  void set(const CGLA::Vec3f& _eye,     // Eye point (camera position)
	   const CGLA::Vec3f& _focus,   // Focus point
	   const CGLA::Vec3f& _up,      // Up vector
	   float fd)                    // Focal distance
  {
#ifndef M_1_PI
    const double M_1_PI = 0.318309886184;
#endif

	eye = _eye;
    focus = _focus;
    up = _up;
    focal_dist = fd;

    line_of_sight = focus - eye;

    // Calculate view plane normal and basis
    vp_normal = normalize(line_of_sight);
    assert(dot(vp_normal, normalize(up)) < 0.99999999);
    vp_axes[0] = normalize(cross(vp_normal, up));
    vp_axes[1] = normalize(cross(vp_axes[0], vp_normal));

    // Account for focal distance
    vp_normal *= focal_dist;
  
    // Calculate field of view (using the pinhole camera model)
    float tmp = atan(1.0/(2.0 * focal_dist));
    fov = 360.0 * M_1_PI * tmp;
  }

  /// Get direction of viewing ray from image coords.
  CGLA::Vec3f get_ray_dir(const CGLA::Vec2f& coords) const
  {
    // vp_normal is multiplied by focal_dist in the constructor
    return vp_normal + vp_axes[0]*coords[0] + vp_axes[1]*coords[1];
  }

  /// Return position of camera.
  const CGLA::Vec3f& get_position() const { return eye; }

  /// Return the ray corresponding to a set of image coords
  Geometry::Ray get_ray(const CGLA::Vec2f& coords) const
  {
    return Geometry::Ray(eye, normalize(get_ray_dir(coords))); 
  }

  float get_fov() const { return fov; }

  float get_focal_dist() const { return focal_dist; }

  // OpenGL

  void glSetPerspective(float width, float height) const
  {
    GLdouble aspect = width/height;    

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov, aspect, focal_dist*NEAR_PLANE, FAR_PLANE);

    glMatrixMode(GL_MODELVIEW);
  }

  void glSetCamera() const
  {
    gluLookAt(eye[0],   eye[1],   eye[2], 
	      focus[0], focus[1], focus[2], 
	      up[0],    up[1],    up[2]);
  }

private:

  CGLA::Vec3f eye, focus, up;
  float focal_dist;
  float fov;
  float phot_rad;

  CGLA::Vec3f line_of_sight;

  // Basis of camera coordinate system (vp - view-plane)
  CGLA::Vec3f vp_normal;
  CGLA::Vec3f vp_axes[2]; 
};

#endif
