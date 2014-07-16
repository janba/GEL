/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

// Created by Bent Dalgaard Larsen, Nov, 2003

#include "Ray.h"
#include "BBox.h"

using namespace std;
using namespace CGLA;

#define my_min(x,y) (x<y?x:y)
#define my_max(x,y) (x>y?x:y)

namespace Geometry
{
  void BBox::intersect_min_max(Ray &ray, double &t_min, double &t_max) const 
  {
    double tx1 = (min_corner[0]-ray.origin[0])/ray.direction[0];
    double ty1 = (min_corner[1]-ray.origin[1])/ray.direction[1];
    double tz1 = (min_corner[2]-ray.origin[2])/ray.direction[2];

    double tx2 = (max_corner[0]-ray.origin[0])/ray.direction[0];
    double ty2 = (max_corner[1]-ray.origin[1])/ray.direction[1];
    double tz2 = (max_corner[2]-ray.origin[2])/ray.direction[2];

    t_min = my_max(my_min(tx1, tx2), my_max(my_min(ty1, ty2), my_min(tz1,tz2)));
    t_max = my_min(my_max(tx1, tx2), my_min(my_max(ty1, ty2), my_max(tz1,tz2)));
  }

  bool BBox::intersect(Ray &ray) 
  {
    double t_min, t_max;
    intersect_min_max(ray, t_min, t_max);
    
    if(t_min <= t_max && t_min < ray.dist && t_min > f_eps) 
      return true;
    else
      return false;
  }

  bool BBox::ray_triangle(Vec3f &ray_start, Vec3f &ray_end, ISectTri &tri) 
  {
    Vec3f origin = ray_start;
    Vec3f direction = ray_end - ray_start;
    double dist = direction.length();
    direction.normalize();
    
    Vec3f p;
    Vec3f q;
    Vec3f s;
    double a, f, u, v, t;
    
    p = cross(direction, Vec3f(tri.edge1));
    a = dot(Vec3f(tri.edge0),p);
    if (a>-f_eps && a<f_eps)
      return false;
    f = 1/a;
    s = origin - Vec3f(tri.point0);
    u = f*dot(s,p);
    if (u<0.0 || u>1.0)
      return false;
    q = cross(s, Vec3f(tri.edge0));
    v = f * dot(direction, q);  
    if (v<0.0 || u+v>1.0)
      return false;
    t = f*dot(Vec3f(tri.edge1), q);
    if (t<0)
      return false;
//	if (t<eps)
//		return false;
    if (t>dist)
      return false;

    return true;
  }


  bool BBox::intersect_edge_box(Vec3f &ray_start, Vec3f &ray_end) 
  {
    Ray test_ray;
    test_ray.origin = ray_start;
    test_ray.direction = ray_end - ray_start;
    test_ray.dist = test_ray.direction.length();
    test_ray.direction.normalize();
    return intersect(test_ray);
  }

  bool BBox::in_interval(double min_limit, double test_value, double max_limit) 
  {
    if (min_limit<=test_value && test_value<=max_limit)
      return true;
    return false;
  }
/*

bool BBox::intersect_triangle_left(ISectTri &tri, double plane, int axis) {
	if (tri.point0[axis]<=plane)
		return true;
	if (tri.point1[axis]<=plane)
		return true;
	if (tri.point2[axis]<=plane)
		return true;
	return false;
}

bool BBox::intersect_triangle_right(ISectTri &tri, double plane, int axis) {
	if (tri.point0[axis]>=plane)
		return true;
	if (tri.point1[axis]>=plane)
		return true;
	if (tri.point2[axis]>=plane)
		return true;
	return false;
}
*/
  bool BBox::intersect_triangle(ISectTri &tri) 
  {
    Vec3f tmin_corner = min_corner - Vec3f(f_eps);
    Vec3f tmax_corner = max_corner + Vec3f(f_eps);

    // Vertex in box test:
    // If any of the triangle vertices are inside the box then 
    // the triangle intersects the box
    if (in_interval(tmin_corner[0],tri.point0[0],tmax_corner[0]) && in_interval(tmin_corner[1],tri.point0[1],tmax_corner[1]) && in_interval(tmin_corner[2],tri.point0[2],tmax_corner[2]))
      return true;
    if (in_interval(tmin_corner[0],tri.point1[0],tmax_corner[0]) && in_interval(tmin_corner[1],tri.point1[1],tmax_corner[1]) && in_interval(tmin_corner[2],tri.point1[2],tmax_corner[2]))
      return true;
    if (in_interval(tmin_corner[0],tri.point2[0],tmax_corner[0]) && in_interval(tmin_corner[1],tri.point2[1],tmax_corner[1]) && in_interval(tmin_corner[2],tri.point2[2],tmax_corner[2]))
      return true;

    // Triangle outside box test:
    // If all of the triangle vertices are outside one of the planes 
    // defining the sides of the box then the triangle can be trivially
    // rejected as outside
    int i;
    for(i=0;i<3;i++)
      if (tri.point0[i]<tmin_corner[i] && tri.point1[i]<tmin_corner[i] && tri.point2[i]<tmin_corner[i])
	return false;
    
    for(i=0;i<3;i++)
      if (tri.point0[i]>tmax_corner[i] && tri.point1[i]>tmax_corner[i] && tri.point2[i]>tmax_corner[i])
	return false;

    // Triangle edges - box intersection test
    if (intersect_edge_box(tri.point0, tri.point1))
      return true;
		
    if (intersect_edge_box(tri.point1, tri.point2))
      return true;

    if (intersect_edge_box(tri.point2, tri.point0))
      return true;

    // Box diagonal - triangle intersection test, 4 tests in total
    Vec3f corner0;
    Vec3f corner1;

    Vec3f tmin_corner_e = tmin_corner;
    Vec3f tmax_corner_e = tmax_corner;

    corner0.set(tmin_corner_e[0],tmin_corner_e[1],tmin_corner_e[2]);
    corner1.set(tmax_corner_e[0],tmax_corner_e[1],tmax_corner_e[2]);
    if (ray_triangle(corner0, corner1, tri))
      return true;

    corner0.set(tmax_corner_e[0],tmin_corner_e[1],tmin_corner_e[2]);
    corner1.set(tmin_corner_e[0],tmax_corner_e[1],tmax_corner_e[2]);
    if (ray_triangle(corner0, corner1, tri))
      return true;

    corner0.set(tmin_corner_e[0],tmax_corner_e[1],tmin_corner_e[2]);
    corner1.set(tmax_corner_e[0],tmin_corner_e[1],tmax_corner_e[2]);
    if (ray_triangle(corner0, corner1, tri))
      return true;

    corner0.set(tmin_corner_e[0],tmin_corner_e[1],tmax_corner_e[2]);
    corner1.set(tmax_corner_e[0],tmax_corner_e[1],tmin_corner_e[2]);
    if (ray_triangle(corner0, corner1, tri))
      return true;

    // None succeded 
    return false;
  }
  
  void BBox::compute_bbox(vector<ISectTri> &isectmesh) 
  {
    min_corner.set(CGLA::BIG,CGLA::BIG,CGLA::BIG);
    max_corner.set(-CGLA::BIG,-CGLA::BIG,-CGLA::BIG);

    for(unsigned int i=0;i<isectmesh.size(); ++i) 
    {
      const ISectTri& tri = isectmesh[i];
      for(int j=0;j<3;j++) {
					if (min_corner[j]>tri.point0[j])
							min_corner[j]=tri.point0[j];
					if (min_corner[j]>tri.point1[j])
							min_corner[j]=tri.point1[j];
					if (min_corner[j]>tri.point2[j])
							min_corner[j]=tri.point2[j];
					if (max_corner[j]<tri.point0[j])
							max_corner[j]=tri.point0[j];
					if (max_corner[j]<tri.point1[j])
							max_corner[j]=tri.point1[j];
					if (max_corner[j]<tri.point2[j])
							max_corner[j]=tri.point2[j];
      }
    }
  }

  double BBox::area() 
  {
    Vec3f size = max_corner - min_corner;
    return size[0]*size[1]*2 + 
      size[1]*size[2]*2 + 
      size[0]*size[2]*2; 
  }
}

