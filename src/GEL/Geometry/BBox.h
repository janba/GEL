/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BBox.h
 * @brief Bounding box class (ancestor of AABox and OBox)
 */

#ifndef __GEOMETRY_BBOX_H__
#define __GEOMETRY_BBOX_H__

#include "../CGLA/Vec3f.h"

#include "Ray.h"

namespace Geometry 
{
  struct ISectTri 
  {
    CGLA::Vec3f point0;
    CGLA::Vec3f point1;
    CGLA::Vec3f point2;
    
    CGLA::Vec3f edge0; // Optimization
    CGLA::Vec3f edge1; // Optimization
    unsigned int mesh_id;
    unsigned int tri_id; // pad to 48 bytes for cache alignment purposes
  };

  struct TriAccel
  {
    // first 16 byte half cache line
    // plane:
    double n_u;  //!< == normal.u / normal.k
    double n_v;  //!< == normal.v / normal.k
    double n_d;  //!< constant of plane equation
    int k;       // projection dimension

    // second 16 byte half cache line
    // line equation for line ac
    double b_nu;
    double b_nv;
    double b_d;
    unsigned int mesh_id;

    // third 16 byte half cache line
    // line equation for line ab
    double c_nu;
    double c_nv;
    double c_d;

    unsigned int tri_id; // pad to 48 bytes for cache alignment purposes
  };

  struct BBox 
  {
    CGLA::Vec3f min_corner;
    CGLA::Vec3f max_corner;

    void intersect_min_max(Ray &ray, double &t_min, double &t_max) const ;
    bool intersect(Ray &ray);
    bool ray_triangle(CGLA::Vec3f &ray_start, CGLA::Vec3f &ray_end, ISectTri &tri);
    bool intersect_edge_box(CGLA::Vec3f &ray_start, CGLA::Vec3f &ray_end);
    bool intersect_triangle(ISectTri &tri);
    bool in_interval(double min_limit, double test_value, double max_limit);
    void compute_bbox(std::vector<ISectTri> &isectmesh);
    bool intersect_triangle_left(ISectTri &tri, double plane, int axis);
    bool intersect_triangle_right(ISectTri &tri, double plane, int axis);
    double area();
  };
}

#endif
