/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Ray.h
 * @brief A ray class for ray tracing.
 */

#ifndef __GEOMETRY_RAY_H__
#define __GEOMETRY_RAY_H__

#include "../CGLA/Vec3i.h"
#include "../CGLA/Vec3f.h"
#include "TriMesh.h"
#include "Material.h"

namespace Geometry 
{
    const double d_eps = 1.0e-12;
    const float f_eps = 1.0e-6f;

    /// Represents a ray as used for ray tracing.
    struct Ray 
    {
        // Constructor
        Ray() 
          : origin(0.0f), direction(0.0f), hit_pos(0.0f), hit_normal(0.0f), 
            has_hit(false), inside(false), did_hit_diffuse(false), dist(CGLA::BIG),
            ior(1.0f), u(0.0f), v(0.0f), trace_depth(0), hit_object(0)
        { }

        Ray(const CGLA::Vec3f& _origin, const CGLA::Vec3f& _direction) 
          : origin(_origin), direction(_direction), hit_pos(0.0f), hit_normal(0.0f), 
            has_hit(false), inside(false), did_hit_diffuse(false), dist(CGLA::BIG),
            ior(1.0f), u(0.0f), v(0.0f), trace_depth(0), hit_object(0)
        { }

        CGLA::Vec3f origin;
        CGLA::Vec3f direction;
        CGLA::Vec3f hit_pos;
        CGLA::Vec3f hit_normal;

        bool has_hit; // Did the ray hit an object
        bool inside;  // Is the ray inside an object
        bool did_hit_diffuse;

        double dist;  // Distance from origin to current intersection
        float ior;    // Current index of refraction for media
        float u, v;   // uv-coordinates on current surface

        int trace_depth;  // Current recursion number
        size_t hit_face_id;
        int id;

        const TriMesh* hit_object;

        const Material* get_hit_material() const
        {
          if(!hit_object)
            return 0;
          return &hit_object->materials[hit_object->mat_idx[hit_face_id]];
        }

        void reset()
        {
            has_hit=false;
            inside = false;
            did_hit_diffuse = false;
            dist = CGLA::BIG;
            ior = 1.0f;
            u=0.0f;
            v=0.0f;
            trace_depth = 0;
            hit_object = 0;
        }

        void compute_position()
        {
            hit_pos = origin + dist*direction;      
        }

        void compute_normal()
        {
            const CGLA::Vec3i& face = hit_object->normals.face(hit_face_id);
            const CGLA::Vec3f& normal0 = hit_object->normals.vertex(face[0]);
            const CGLA::Vec3f& normal1 = hit_object->normals.vertex(face[1]);
            const CGLA::Vec3f& normal2 = hit_object->normals.vertex(face[2]);
            hit_normal = normalize(normal0*(1 - u - v) + normal1*u + normal2*v);      
        }

        void reflect(const CGLA::Vec3f &normal)
        {
            assert(dot(direction, normal) < 0.0f);
            direction = normal*2.0f*dot(-direction,normal) + direction;      
        }

        void refract(const CGLA::Vec3f& normal, float new_ior)
        {
            float ref_ratio = ior/new_ior;
            float cos_N_I = dot(normal, direction);
            CGLA::Vec3f norm(normal);

            if(cos_N_I > 0.0f)
            {
                norm = -norm;
                cos_N_I = dot(norm, direction);
            }

            float selector = 1+(ref_ratio*ref_ratio)*(cos_N_I*cos_N_I - 1);

            if(selector > 0.0f) 
            {
                direction = norm*(ref_ratio*(-cos_N_I) - sqrt(selector)) 
                            + direction*ref_ratio;
                ior = new_ior;
            } 
            else 
                // Total internal reflection.
                reflect(normal);
        }

        bool cond_set_parameter(float t, float _u, float _v, const Geometry::TriMesh* mesh, size_t idx)
        {
            if(t < dist)
            {
                dist = t;
                u = _u;
                v = _v;
                hit_object = mesh;
                hit_face_id = idx;
                has_hit = true;
                return true;
            }
            return false;
        }
    };
}
#endif

